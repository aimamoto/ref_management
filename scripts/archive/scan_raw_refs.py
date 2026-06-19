import re
import csv
import os
import sys
import argparse
import time
from docx import Document
from Bio import Entrez
from rapidfuzz import fuzz

# --- CONFIGURATION ---
Entrez.email = os.environ.get("NCBI_EMAIL", "aimamoto@uchicago.edu")
Entrez.api_key = os.environ.get("NCBI_API_KEY", None)

if Entrez.api_key is None:
    print("Tip: Set NCBI_API_KEY environment variable to avoid request limits.")

# --- HELPERS ---
def clean_word(word):
    return re.sub(r'[^\w]', '', word)

def resolve_doi_to_pmid(doi):
    try:
        clean = doi.rstrip('.').strip()
        handle = Entrez.esearch(db="pubmed", term=f"{clean}[DOI]", retmax=1)
        r = Entrez.read(handle)
        handle.close()
        return r['IdList'][0] if r['IdList'] else None
    except: return None

def extract_citation_parts(text):
    """Extracts Year, Author, and Title Snippet for search."""
    data = {"year": "", "author": "", "title_snippet": ""}
    
    # 1. Year
    year_match = re.search(r'\b(19|20)\d{2}\b', text)
    if year_match:
        data["year"] = year_match.group(0)
    
    # 2. Author
    words = text.split()
    for w in words:
        clean = clean_word(w)
        if clean.isdigit() or len(clean) < 2: continue
        if clean.lower() in ['et', 'al', 'in', 'the', 'pmid', 'doi', 'vol', 'no']: continue
        data["author"] = clean
        break
    
    # 3. Title Snippet
    if data["year"]:
        pattern = re.escape(data["year"]) + r'\)?\.?\s+([^.]+)\.'
        title_match = re.search(pattern, text)
        if title_match:
            data["title_snippet"] = title_match.group(1).strip()
        else:
            start = text.find(data["year"]) + 4
            data["title_snippet"] = text[start:start+100]
    else:
        data["title_snippet"] = text[:80]
            
    return data

def search_pubmed_by_metadata(parts):
    candidates = []
    # Strategy A: Title Search
    if len(parts["title_snippet"]) > 15:
        clean_title = re.sub(r'[^\w\s]', '', parts["title_snippet"])
        short_title = " ".join(clean_title.split()[:8])
        query = f"{short_title}[Title]"
        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=3)
            r = Entrez.read(handle)
            handle.close()
            candidates.extend(r['IdList'])
        except: pass
    
    # Strategy B: Author + Year
    if not candidates and parts["author"] and parts["year"]:
        query = f"{parts['author']}[1au] AND {parts['year']}[pdat]"
        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=5)
            r = Entrez.read(handle)
            handle.close()
            candidates.extend(r['IdList'])
        except: pass
        
    return list(set(candidates))

def parse_record(record):
    title = re.sub(r'<[^<]+?>', '', record.get('Title', ''))
    authors = record.get('AuthorList', [])
    first_author = authors[0] if authors else "Unknown"
    pub_date = record.get('PubDate', '')
    year_match = re.search(r'\d{4}', pub_date)
    year = year_match.group(0) if year_match else "Unknown"
    
    doi = ""
    for aid in record.get('ArticleIds', {}).items():
        if aid[0] == 'doi': doi = aid[1]
    if not doi:
        ids = record.get('ArticleIds', {})
        if 'doi' in ids: doi = ids['doi']

    is_retracted = "Retracted Publication" in record.get('PubTypeList', [])

    return {
        "id": record.get('Id'),
        "first_author": first_author,
        "year": year,
        "title": title,
        "doi": doi,
        "is_retracted": is_retracted
    }

def batch_fetch_pubmed(pmid_list):
    if not pmid_list: return {}
    BATCH_SIZE = 200
    fetched_data = {}
    unique_pmids = list(set(pmid_list))
    
    print(f"\nFetching metadata for {len(unique_pmids)} unique papers...")
    for i in range(0, len(unique_pmids), BATCH_SIZE):
        chunk = unique_pmids[i : i + BATCH_SIZE]
        ids_str = ",".join(chunk)
        try:
            handle = Entrez.esummary(db="pubmed", id=ids_str, retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            for record in records:
                data = parse_record(record)
                fetched_data[data['id']] = data
            time.sleep(0.5) 
        except Exception as e:
            print(f"Warning: Batch fetch error: {e}")
            
    return fetched_data

# --- CORE LOGIC ---
def analyze_document(file_path):
    print(f"\nReading document: {file_path}...")
    
    all_paragraphs = []
    ext = os.path.splitext(file_path)[1].lower()
    
    try:
        if ext == '.docx':
            doc = Document(file_path)
            all_paragraphs = [p.text.strip() for p in doc.paragraphs if p.text.strip()]
        elif ext == '.txt':
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                lines = f.readlines()
                current_block = []
                for line in lines:
                    if line.strip():
                        # Detect list items or new paragraphs
                        if re.match(r'^\[?\d+[\.\)\]]', line.strip()):
                             all_paragraphs.append(line.strip())
                        else:
                             current_block.append(line.strip())
                    else:
                        if current_block:
                            all_paragraphs.append(" ".join(current_block))
                            current_block = []
                if current_block: all_paragraphs.append(" ".join(current_block))
        else:
            print(f"ERROR: Unsupported file format {ext}")
            return [], []
    except Exception as e:
        print(f"ERROR: Could not open file. {e}")
        return [], []

    # --- SECTION DETECTION ---
    header_regex = re.compile(r'^\s*(?:[0-9]+\.?\s*)?(?:REFERENCES|BIBLIOGRAPHY|LITERATURE CITED)\s*$', re.IGNORECASE)
    
    start_index = 0
    header_found = False
    for i, p in enumerate(all_paragraphs):
        if header_regex.match(p):
            start_index = i + 1
            header_found = True
            print(f" -> Found Reference Section header at paragraph {i+1}. Processing subsequent text.")
            break
    
    if not header_found:
        print(" -> No 'References' header found. Assuming file is just a list of references.")

    paragraphs_to_check = all_paragraphs[start_index:]
    
    pmid_pattern = re.compile(r'PMID:?\s*(\d+)', re.IGNORECASE)
    doi_pattern = re.compile(r'\b(10\.\d{4,9}/[-._;()/:a-zA-Z0-9]+)\b')

    items_to_process = []
    ids_to_fetch_metadata = set()

    print(f"Scanning {len(paragraphs_to_check)} reference candidates for IDs...")
    
    for text in paragraphs_to_check:
        if len(text) < 15: continue 

        found_pmids = [m.group(1) for m in pmid_pattern.finditer(text)]
        found_dois = [m.group(1).rstrip('.') for m in doi_pattern.finditer(text)]
        
        entry = {
            "original_text": text,
            "found_pmid": found_pmids[0] if found_pmids else None,
            "found_doi": found_dois[0] if found_dois else None,
            "target_pmid": None, 
            "status": "PENDING",
            "action_log": []
        }

        # SCENARIO 1: Has PMID
        if entry["found_pmid"]:
            entry["target_pmid"] = entry["found_pmid"]
            ids_to_fetch_metadata.add(entry["target_pmid"])
            
        # SCENARIO 2: No PMID, but Has DOI
        elif entry["found_doi"]:
            # Note: We resolve one by one here to get the ID for the batch fetch later
            resolved_pmid = resolve_doi_to_pmid(entry["found_doi"])
            if resolved_pmid:
                entry["target_pmid"] = resolved_pmid
                entry["action_log"].append("Fetched PMID via DOI")
                ids_to_fetch_metadata.add(resolved_pmid)
            else:
                entry["status"] = "FAIL_DOI_LOOKUP"

        # SCENARIO 3: Missing Both
        else:
            if re.search(r'\b(19|20)\d{2}\b', text):
                parts = extract_citation_parts(text)
                candidates = search_pubmed_by_metadata(parts)
                
                best_match = None
                best_score = 0
                if candidates:
                    cand_meta = batch_fetch_pubmed(candidates)
                    for c_id, c_data in cand_meta.items():
                        score = fuzz.token_set_ratio(c_data['title'], parts['title_snippet'])
                        if score > 85 and score > best_score:
                            best_score = score
                            best_match = c_id
                
                if best_match:
                    entry["target_pmid"] = best_match
                    entry["action_log"].append("Found via Search")
                    ids_to_fetch_metadata.add(best_match)
        
        items_to_process.append(entry)

    pubmed_db = batch_fetch_pubmed(list(ids_to_fetch_metadata))
    
    csv_results = []
    txt_lines = []
    
    # --- CHATTY OUTPUT ---
    print("\n" + "="*85)
    print(f"{'ID / Ref Start':<18} | {'Status':<18} | {'Action / Note'}")
    print("-" * 85)

    for item in items_to_process:
        text = item['original_text']
        pmid = item['target_pmid']
        meta = pubmed_db.get(pmid)
        
        status = "OK"
        notes = item['action_log']
        display_id = pmid if pmid else "(No ID)"
        
        if not pmid:
            if item['status'] == "FAIL_DOI_LOOKUP":
                status = "MANUAL_CHECK"
                notes.append("DOI invalid/not in PubMed")
            else:
                if re.search(r'\b(19|20)\d{2}\b', text):
                    status = "NOT_FOUND"
                    display_id = "UNKNOWN"
                else:
                    status = "IGNORED" 
        
        elif meta:
            parts = extract_citation_parts(text)
            if meta['is_retracted']:
                status = "!! RETRACTED !!"
            
            elif "Found via Search" not in str(notes):
                # Verify match
                title_score = fuzz.token_set_ratio(meta['title'], parts['title_snippet'])
                if title_score < 70:
                    status = "MISMATCH"
                    notes.append("ID mismatch with text")

        # --- OUTPUT GENERATION ---
        if status == "IGNORED": continue 

        clean_line = text.replace('\n', ' ').strip()
        
        if meta:
            # 1. Inject Missing PMID
            if not item['found_pmid']:
                clean_line += f" [PMID: {pmid}]"
                if "Fetched PMID" not in str(notes) and "Found via Search" not in str(notes):
                     notes.append("Added Missing PMID")

            # 2. Inject Missing DOI
            if meta['doi'] and meta['doi'] not in clean_line:
                clean_line += f" [doi: {meta['doi']}]"
                notes.append("Added DOI")
            
            if "Found via Search" in str(notes):
                 clean_line += f" [Verified: {meta['first_author']}, {meta['year']}]"

        # 4. Warnings
        if status == "!! RETRACTED !!":
            clean_line += " [!! RETRACTED PUBLICATION !!]"
        elif status == "MISMATCH":
            clean_line += " [!! ID MISMATCH - CHECK MANUALLY !!]"
        elif status == "NOT_FOUND":
            clean_line += " [!! NO IDENTIFIER FOUND !!]"

        # Print formatted row
        action_str = ", ".join(notes) if notes else "Verified"
        if len(display_id) > 15: display_id = display_id[:15]
        print(f"{display_id:<18} | {status:<18} | {action_str}")

        csv_results.append({
            "Original_Text": text[:60] + "...",
            "Final_PMID": pmid if pmid else "N/A",
            "Status": status,
            "Notes": "; ".join(notes)
        })
        txt_lines.append(clean_line)
    
    print("="*85)
    return csv_results, txt_lines

def save_outputs(csv_data, txt_lines, base_name):
    if not csv_data: return
    csv_name = f"{base_name}_verification_report.csv"
    txt_name = f"{base_name}_corrected_refs.txt"

    with open(csv_name, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=csv_data[0].keys())
        writer.writeheader()
        writer.writerows(csv_data)
    
    with open(txt_name, 'w', encoding='utf-8') as f:
        f.write(f"--- Processed References for {base_name} ---\n\n")
        for line in txt_lines:
            f.write(line + "\n\n")
    print(f"\nSuccess! Saved to {csv_name} and {txt_name}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Auto-verify references. Detects 'References' section.")
    parser.add_argument("file", help="Path to the .docx or .txt file", nargs='?', default=None)
    args = parser.parse_args()
    fname = args.file

    if not fname:
        fname = input("Enter filename (.docx or .txt): ").strip().strip('"').strip("'")
    
    if not os.path.exists(fname):
        print(f"Error: File '{fname}' not found.")
        sys.exit(1)

    base = re.split(r'[ _]', os.path.splitext(os.path.basename(fname))[0])[0] or "scan"
    c, t = analyze_document(fname)
    save_outputs(c, t, base)
