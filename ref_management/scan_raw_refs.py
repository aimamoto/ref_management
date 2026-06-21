import re
import csv
import os
import sys
import argparse
import time
import requests
from docx import Document
from Bio import Entrez
from rapidfuzz import fuzz

# --- CONFIGURATION ---
Entrez.email = os.environ.get("NCBI_EMAIL", None)
Entrez.api_key = os.environ.get("NCBI_API_KEY", None)

if Entrez.api_key is None:
    print("Tip: Set NCBI_API_KEY environment variable to avoid request limits.", file=sys.stderr)

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
    except Exception:
        return None

def resolve_doi_crossref(doi):
    """Fallback for DOIs not indexed in PubMed (e.g. Stats journals, old issues)."""
    print(f"      ...PubMed missed {doi}. Checking Crossref Master Database...", end='\r')
    try:
        clean = doi.rstrip('.').strip()
        url = f"https://api.crossref.org/works/{clean}"
        resp = requests.get(url, timeout=5)
        if resp.status_code == 200:
            data = resp.json().get('message', {})
            title = data.get('title', [''])[0]
            author_list = data.get('author', [{}])
            first_author = author_list[0].get('family', 'Unknown') if author_list else "Unknown"
            
            year = "Unknown"
            pub = data.get('published-print') or data.get('published-online')
            if pub and 'date-parts' in pub:
                year = str(pub['date-parts'][0][0])
            
            return {
                "id": f"CR_{clean}",
                "first_author": first_author,
                "year": year,
                "title": title,
                "doi": clean,
                "is_retracted": False
            }
    except Exception:
        pass
    return None

def extract_citation_parts(text):
    data = {"year": "", "author": "", "title_snippet": ""}
    year_iter = list(re.finditer(r'\b(19|20)\d{2}[a-z]?\b', text))
    
    if year_iter:
        selected_year_match = year_iter[0]
        for m in year_iter:
            if m.start() > 0 and text[m.start()-1] == '(':
                selected_year_match = m
                break
                
        raw_year = selected_year_match.group(0)
        data["year"] = re.sub(r'[a-z]', '', raw_year)
        
        # Determine if Author-Year or Vancouver style based on Year placement
        if selected_year_match.start() > len(text) / 2:
            snippet = text[:selected_year_match.start()].strip()
        else:
            snippet = text[selected_year_match.end():].strip()
            
        # URL/DOI Scrubber (Prevents them from leaking into the Title)
        snippet = re.sub(r'(?i)https?://\S+', '', snippet)
        snippet = re.sub(r'(?i)doi:?\s*10\.\d{4,9}/[-._;()/:a-zA-Z0-9<>\[\]]+', '', snippet)
        data["title_snippet"] = re.sub(r'^[\s\)\.,:;-]+|[\s\.,:;-]+$', '', snippet)
        
        # Improved Raw Author Parsing
        author_raw = text[:selected_year_match.start()].split('.')[0].strip()
        author_raw = re.sub(r'^\[?\d+\]?\s*', '', author_raw) 
        if len(author_raw) > 3:
            data["author"] = author_raw
    else:
        snippet = text
        snippet = re.sub(r'(?i)https?://\S+', '', snippet)
        snippet = re.sub(r'(?i)doi:?\s*10\.\d{4,9}/[-._;()/:a-zA-Z0-9<>\[\]]+', '', snippet)
        data["title_snippet"] = snippet.strip()

    if not data["author"]:
        skip_words = {'et', 'al', 'in', 'the', 'pmid', 'doi', 'vol', 'no', 'and', '&', 'eds', 'editor', 'page', 'pp', 'references'}
        for w in text.split():
            clean = clean_word(w)
            if len(clean) < 2 or clean.isdigit(): continue
            if clean.lower() in skip_words: continue
            if clean[0].isupper():
                data["author"] = clean
                break
            
    return data

def search_pubmed_by_metadata(parts):
    candidates = []
    if len(parts["title_snippet"]) > 10:
        clean_title = re.sub(r'[^\w\s]', '', parts["title_snippet"])
        short_title = " ".join(clean_title.split()[:8])
        try:
            handle = Entrez.esearch(db="pubmed", term=f"{short_title}[Title]", retmax=3)
            r = Entrez.read(handle)
            handle.close()
            if r['IdList']: candidates.extend(r['IdList'])
        except Exception: pass
    
    if not candidates and parts["author"] and parts["year"]:
        try:
            handle = Entrez.esearch(db="pubmed", term=f"{parts['author'].split()[0]}[1au] AND {parts['year']}[pdat]", retmax=5)
            r = Entrez.read(handle)
            handle.close()
            if r['IdList']: candidates.extend(r['IdList'])
        except Exception: pass
        
    return list(set(candidates))

def parse_record(record):
    title = re.sub(r'<[^<]+?>', '', record.get('Title', ''))
    authors = record.get('AuthorList', [])
    first_author = authors[0] if authors else "Unknown"
    pub_date = record.get('PubDate', '')
    year_match = re.search(r'\d{4}', pub_date)
    
    doi = ""
    for aid in record.get('ArticleIds', {}).items():
        if aid[0] == 'doi': doi = aid[1]
    if not doi and 'doi' in record.get('ArticleIds', {}):
        doi = record['ArticleIds']['doi']

    return {
        "id": str(record.get('Id')),
        "first_author": first_author,
        "year": year_match.group(0) if year_match else "Unknown",
        "title": title,
        "doi": doi,
        "is_retracted": "Retracted Publication" in record.get('PubTypeList', [])
    }

def batch_fetch_pubmed(pmid_list):
    if not pmid_list: return {}
    fetched_data = {}
    unique_pmids = list(set(pmid_list))
    
    print(f"\nFetching metadata for {len(unique_pmids)} unique PubMed papers...")
    for i in range(0, len(unique_pmids), 200):
        chunk = unique_pmids[i : i + 200]
        try:
            handle = Entrez.esummary(db="pubmed", id=",".join(chunk), retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            for record in records:
                data = parse_record(record)
                fetched_data[data['id']] = data
            time.sleep(0.5) 
        except Exception as e:
            pass
            
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
                        if re.match(r'^\[?\d+[\.\)\]]', line.strip()):
                             if current_block: all_paragraphs.append(" ".join(current_block))
                             current_block = [line.strip()]
                             all_paragraphs.append(line.strip())
                             current_block = []
                        else:
                             current_block.append(line.strip())
                    else:
                        if current_block:
                            all_paragraphs.append(" ".join(current_block))
                            current_block = []
                if current_block: all_paragraphs.append(" ".join(current_block))
        else:
            print(f"ERROR: Unsupported format {ext}"); return [], [], []
    except Exception as e:
        print(f"ERROR: Could not open file. {e}"); return [], [], []

    header_regex = re.compile(r'^\s*(?:[0-9]+\.?\s*)?(?:REFERENCES|BIBLIOGRAPHY|LITERATURE CITED|WORKS CITED)\s*$', re.IGNORECASE)
    start_index, header_found = 0, False
    for i, p in enumerate(all_paragraphs):
        if header_regex.match(p):
            start_index, header_found = i + 1, True
            print(f" -> Found Reference Section header. Processing subsequent text.")
            break
    if not header_found: print(" -> No 'References' header found. Assuming file is just a list.")

    paragraphs_to_check = all_paragraphs[start_index:]
    pmid_pattern = re.compile(r'PMID:?\s*(\d+)', re.IGNORECASE)
    
    # NEW DOI PATTERN: Supports <, >, and brackets used in Wiley/SICI DOIs
    doi_pattern = re.compile(r'\b(10\.\d{4,9}/[-._;()/:a-zA-Z0-9<>\[\]]+)')

    items_to_process = []
    ids_to_fetch_metadata = set()

    print(f"Scanning {len(paragraphs_to_check)} reference candidates...")
    
    for idx, text in enumerate(paragraphs_to_check):
        if len(text) < 15: continue 

        num_match = re.match(r'^\s*\[?(\d+)[\.\)\]]', text)
        ref_num = num_match.group(1) if num_match else str(idx + 1)

        found_pmids = [m.group(1) for m in pmid_pattern.finditer(text)]
        found_dois = [m.group(1).rstrip('.') for m in doi_pattern.finditer(text)]
        
        entry = {
            "ref_num": ref_num,
            "original_text": text,
            "found_pmid": found_pmids[0] if found_pmids else None,
            "found_doi": found_dois[0] if found_dois else None,
            "target_pmid": None,
            "cr_meta": None,
            "status": "PENDING",
            "action_log": []
        }

        if entry["found_pmid"]:
            entry["target_pmid"] = entry["found_pmid"]
            ids_to_fetch_metadata.add(entry["target_pmid"])
        elif entry["found_doi"]:
            resolved_pmid = resolve_doi_to_pmid(entry["found_doi"])
            if resolved_pmid:
                entry["target_pmid"] = resolved_pmid
                entry["action_log"].append("Fetched PMID via DOI")
                ids_to_fetch_metadata.add(resolved_pmid)
            else:
                cr_meta = resolve_doi_crossref(entry["found_doi"])
                if cr_meta:
                    entry["cr_meta"] = cr_meta
                    entry["action_log"].append("Resolved via Crossref")
                else:
                    entry["status"] = "FAIL_DOI_LOOKUP"
        else:
            if re.search(r'\b(19|20)\d{2}[a-z]?\b', text):
                parts = extract_citation_parts(text)
                candidates = search_pubmed_by_metadata(parts)
                best_match, best_score = None, 0
                if candidates:
                    cand_meta = batch_fetch_pubmed(candidates)
                    for c_id, c_data in cand_meta.items():
                        score = fuzz.token_set_ratio(c_data['title'], parts['title_snippet'])
                        if score > 80 and score > best_score:
                            best_score, best_match = score, c_id
                if best_match:
                    entry["target_pmid"] = best_match
                    entry["action_log"].append("Found via Search")
                    ids_to_fetch_metadata.add(best_match)
        items_to_process.append(entry)

    pubmed_db = batch_fetch_pubmed(list(ids_to_fetch_metadata))
    csv_results, txt_lines, bib_entries = [], [], []
    
    print("\n" + "="*85)
    for item in items_to_process:
        text, pmid, cr_meta = item['original_text'], item['target_pmid'], item['cr_meta']
        meta = pubmed_db.get(pmid) if pmid else cr_meta
        status, notes, display_id = "OK", item['action_log'], pmid if pmid else (item.get('found_doi') if cr_meta else "(No ID)")
        
        if not meta:
            if item['status'] == "FAIL_DOI_LOOKUP": 
                status, _ = "MANUAL_CHECK", notes.append("DOI not in PubMed/Crossref")
            else: 
                status, display_id = ("NOT_FOUND", "UNKNOWN") if re.search(r'\b(19|20)\d{2}\b', text) else ("IGNORED", display_id)
        else:
            title_clean = re.sub(r'[^\w\s]', '', meta['title'].lower())
            text_clean = re.sub(r'[^\w\s]', '', text.lower())
            title_score = fuzz.partial_ratio(title_clean, text_clean)
            
            if meta['is_retracted']: 
                status = "!! RETRACTED !!"
            elif "Found via Search" not in str(notes) and "Crossref" not in str(notes):
                if title_score < 65:
                    status, _ = "MISMATCH", notes.append(f"Title Score {title_score}")

        if status == "IGNORED": continue 

        csv_results.append({"Original_Text": text[:60] + "...", "Final_ID": display_id, "Status": status, "Notes": "; ".join(notes)})
        
        bib_lines = [f"@article{{{item['ref_num']},"]
        if pmid: bib_lines.append(f"  pmid = {{{pmid}}},")
        if meta and meta.get('doi'): bib_lines.append(f"  doi = {{{meta['doi']}}},")
        elif item.get('found_doi'): bib_lines.append(f"  doi = {{{item['found_doi']}}},")
        
        if meta:
            bib_lines.extend([f"  title = {{{meta['title']}}},", f"  author = {{{meta['first_author']} et al.}},", f"  year = {{{meta['year']}}}"])
        else:
            parts = extract_citation_parts(text)
            bib_lines.append(f"  title = {{{parts['title_snippet']}}},")
            if parts['author']: bib_lines.append(f"  author = {{{parts['author']}}},")
            if parts['year']: bib_lines.append(f"  year = {{{parts['year']}}}")
        bib_lines.append("}")
        bib_entries.append("\n".join(bib_lines))
        
        print(f"{display_id[:15]:<18} | {status:<18} | {', '.join(notes) if notes else 'Verified'}")
    
    print("="*85)
    return csv_results, txt_lines, bib_entries

def save_outputs(csv_data, txt_lines, bib_entries, base_name):
    if not csv_data: return
    csv_name, bib_name = f"{base_name}_verification_report.csv", f"{base_name}_extracted.bib"

    with open(csv_name, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=csv_data[0].keys())
        writer.writeheader()
        writer.writerows(csv_data)
        
    with open(bib_name, 'w', encoding='utf-8') as f:
        f.write("\n\n".join(bib_entries))
        
    print(f"\nSuccess! Saved outputs:")
    print(f"  [Report] -> {csv_name}")
    print(f"  [BibTeX] -> {bib_name}  <-- USE THIS FOR arm-verify")

def main():
    parser = argparse.ArgumentParser(description="Scan docx, output reports and a mapped .bib file.")
    parser.add_argument("file", help="Path to .docx", nargs='?', default=None)
    args = parser.parse_args()
    fname = args.file

    if not fname: fname = input("Enter filename (.docx): ").strip().strip('"').strip("'")
    if not fname or not os.path.exists(fname): sys.exit(1)

    base = re.split(r'[ _]', os.path.splitext(os.path.basename(fname))[0])[0] or "scan"
    c, t, b = analyze_document(fname)
    save_outputs(c, t, b, base)

if __name__ == "__main__":
    main()
