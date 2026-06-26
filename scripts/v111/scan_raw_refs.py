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

# --- POST-REFERENCE BOUNDARY ---
POST_REF_PATTERN = re.compile(r'^\s*(?:Tables?|Figures?|Figure Legends?|Supplementary.*?|Appendices|Data Availability|Acknowledgements?|Author Contributions?|Funding|Conflict(?:s)? of Interest|Competing Interests?|(?:Table|Figure|Fig\.?)\s*\d+.*)$', re.IGNORECASE)

# --- HELPERS ---
def safe_entrez_call(func, *args, **kwargs):
    """Respects NCBI Entrez rate limits by dynamically throttling requests."""
    if not Entrez.api_key:
        time.sleep(0.35)
    else:
        time.sleep(0.1)
    return func(*args, **kwargs)

def clean_word(word):
    return re.sub(r'[^\w]', '', word)

def clean_for_ratio(text: str) -> str:
    """Normalizes titles by stripping punctuation and separating hyphens/slashes to preserve terms."""
    if not text:
        return ""
    text = str(text).replace('{', '').replace('}', '')
    text = text.replace('-', ' ').replace('–', ' ').replace('/', ' ')
    return re.sub(r'[^\w\s]', '', text.lower()).strip()

def resolve_doi_to_pmid(doi):
    try:
        clean = doi.rstrip('.').strip()
        handle = safe_entrez_call(Entrez.esearch, db="pubmed", term=f"{clean}[DOI]", retmax=1)
        r = Entrez.read(handle)
        handle.close()
        return r['IdList'][0] if r['IdList'] else None
    except Exception:
        return None

def resolve_doi_crossref(doi):
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

def clean_author_field(author_str: str) -> str:
    """Normalizes a raw author list string into a BibTeX-compatible 'and' separated format."""
    if not author_str:
        return "Unknown"
    author_str = re.sub(r'^\[?\d+\]?\s*', '', author_str).strip()
    author_str = re.sub(r'\s+&\s+', ' and ', author_str)
    
    if ' and ' in author_str:
        parts = [p.strip() for p in author_str.split(' and ') if p.strip()]
    else:
        raw_parts = [p.strip() for p in author_str.split(',') if p.strip()]
        parts = []
        is_paired = False
        if len(raw_parts) > 1:
            second_part_clean = re.sub(r'[.\s]', '', raw_parts[1])
            if len(second_part_clean) <= 3 and second_part_clean.isalpha():
                is_paired = True
                
        if is_paired:
            for idx in range(0, len(raw_parts), 2):
                if idx + 1 < len(raw_parts):
                    parts.append(f"{raw_parts[idx]}, {raw_parts[idx+1]}")
                else:
                    parts.append(raw_parts[idx])
        else:
            parts = raw_parts
            
    clean_parts = []
    for p in parts:
        p_clean = p.strip()
        if p_clean:
            clean_parts.append(p_clean)
            
    if not clean_parts:
        return "Unknown"
    return " and ".join(clean_parts)

def get_author_surname(author_str: str) -> str:
    """Intelligently extracts the author's longest word from the first author block."""
    if not author_str:
        return ""
    author_str = re.sub(r'^\[?\d+\]?\s*', '', author_str.strip())
    first_block = author_str.split(',')[0].split(' and ')[0].strip()
    words = [re.sub(r'[^A-Za-z]', '', w) for w in first_block.split()]
    valid_words = [w for w in words if len(w) > 1 and w.lower() not in ['jr', 'sr', 'iii', 'ii']]
    if valid_words:
        return max(valid_words, key=len)
    return first_block[:10].strip()

def verify_match_precision(draft_title: str, fetched_title: str, draft_author: str, fetched_first_author: str) -> bool:
    title_score = fuzz.ratio(clean_for_ratio(draft_title), clean_for_ratio(fetched_title))
    if title_score < 75:
        return False
    if draft_author and fetched_first_author:
        draft_surname = get_author_surname(draft_author)
        fetched_surname = get_author_surname(fetched_first_author)
        if draft_surname and fetched_surname:
            author_score = fuzz.ratio(clean_for_ratio(draft_surname), clean_for_ratio(fetched_surname))
            if author_score < 60:
                return False
    return True

def find_publication_year(text: str):
    candidates = list(re.finditer(r'\b((?:19|20)\d{2})[a-z]?\b', text))
    if not candidates:
        return None
        
    best_candidate = None
    best_score = -9999
    
    for m in candidates:
        year_val = int(m.group(1))
        if not (1850 <= year_val <= 2028):
            continue
            
        score = 0
        start_idx = m.start()
        end_idx = m.end()
        
        has_left_paren = start_idx > 0 and text[start_idx - 1] in '(['
        has_right_paren = end_idx < len(text) and text[end_idx] in ')]'
        if has_left_paren and has_right_paren:
            score += 150
        elif has_left_paren or has_right_paren:
            score += 50
            
        prefix = text[max(0, start_idx-2):start_idx]
        suffix = text[end_idx:end_idx+2]
        if re.search(r'[\s,;:\.]', prefix): score += 10
        if re.search(r'[\s,;:\.]', suffix): score += 10
            
        char_before = text[start_idx-1] if start_idx > 0 else ''
        char_after = text[end_idx] if end_idx < len(text) else ''
        if char_before in '-–' or char_after in '-–': score -= 100
        if re.match(r'^(?:bp|kb|mb|hz|rpm|g|v|c|s|h|d|m|nm|um|mm|cm)\b', suffix, re.IGNORECASE): score -= 80

        rel_pos = start_idx / len(text)
        if rel_pos < 0.40: score += 20
        elif rel_pos > 0.75: score += 30
            
        if score > best_score:
            best_score = score
            best_candidate = m
            
    return best_candidate

def extract_citation_parts(text: str) -> dict:
    data = {"year": "", "author": "", "title_snippet": ""}
    year_match = find_publication_year(text)
    
    if year_match:
        data["year"] = year_match.group(1)
        start_idx = year_match.start()
        end_idx = year_match.end()
        
        rel_pos = start_idx / len(text)
        is_near_front = (rel_pos < 0.4) or (start_idx > 0 and text[start_idx - 1] in '([')
        
        if is_near_front:
            author_raw = text[:start_idx].strip()
            snippet = text[end_idx:].strip()
            snippet = re.sub(r'(?i)https?://\S+', '', snippet)
            snippet = re.sub(r'(?i)doi:?\s*10\.\d{4,9}/[-._;()/:a-zA-Z0-9<>\[\]]+', '', snippet)
            # Correctly isolate title from journal info
            title_candidate = re.split(r'[\.\?\!]\s+(?=[A-Z])', snippet)[0].strip()
        else:
            text_before = text[:start_idx].strip()
            text_before = re.sub(r'(?i)https?://\S+', '', text_before)
            text_before = re.sub(r'(?i)doi:?\s*10\.\d{4,9}/[-._;()/:a-zA-Z0-9<>\[\]]+', '', text_before)
            parts = re.split(r'\.\s+', text_before)
            
            if len(parts) >= 2:
                author_raw = parts[0].strip()
                title_candidate = parts[1].strip()
            else:
                author_raw = text_before
                title_candidate = text_before
                
        author_raw = re.sub(r'[\(\[\s,;:-]+$', '', author_raw)
        author_raw = re.sub(r'^\[?\d+\]?\s*', '', author_raw)
        if len(author_raw) > 3:
            data["author"] = author_raw
            
        title_candidate = re.sub(r'^[\s\)\.,:;-]+|[\s\.,:;-]+$', '', title_candidate)
        data["title_snippet"] = title_candidate
        
    else:
        snippet = text
        snippet = re.sub(r'(?i)https?://\S+', '', snippet)
        snippet = re.sub(r'(?i)doi:?\s*10\.\d{4,9}/[-._;()/:a-zA-Z0-9<>\[\]]+', '', snippet)
        title_candidate = re.split(r'[\.\?\!]\s+(?=[A-Z])', snippet)[0].strip()
        data["title_snippet"] = title_candidate.strip()

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
    surname = get_author_surname(parts["author"])
    
    if len(parts["title_snippet"]) > 10:
        clean_title = parts["title_snippet"].replace('-', ' ').replace('–', ' ').replace('/', ' ')
        clean_title = re.sub(r'[^\w\s]', '', clean_title)
        short_title = " ".join(clean_title.split()[:8])
        
        # Free-text NLP search
        query = short_title
        if surname:
            query += f" {surname}[1au]"
        else:
            query += "[Title]"
            
        try:
            handle = safe_entrez_call(Entrez.esearch, db="pubmed", term=query, retmax=5)
            r = Entrez.read(handle)
            handle.close()
            if r['IdList']: candidates.extend(r['IdList'])
        except Exception: pass
    
    if not candidates and surname and parts["year"]:
        try:
            handle = safe_entrez_call(Entrez.esearch, db="pubmed", term=f"{surname}[1au] AND {parts['year']}[pdat]", retmax=5)
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
    if not doi and record.get('DOI'):
        doi = record['DOI']

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
            handle = safe_entrez_call(Entrez.esummary, db="pubmed", id=",".join(chunk), retmode="xml")
            records = Entrez.read(handle)
            handle.close()
            for record in records:
                data = parse_record(record)
                fetched_data[data['id']] = data
        except Exception:
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
            for p_node in doc.element.body.xpath('.//w:p'):
                text = "".join(node.text for node in p_node.xpath('.//w:t') if node.text).strip()
                if text: 
                    all_paragraphs.append(text)
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

    header_regex = re.compile(r'^\s*(?:[0-9]+\.?\s*)?(?:REFERENCES?|REFERENCE LIST|BIBLIOGRAPHY|LITERATURE CITED|WORKS CITED)\s*$', re.IGNORECASE)
    start_index, header_found = 0, False
    for i, p in enumerate(all_paragraphs):
        if header_regex.match(p):
            start_index, header_found = i + 1, True
            print(f" -> Found Reference Section header. Processing subsequent text.")
            break
            
    if not header_found: 
        print(" -> No 'References' header found. Scanning entire document from the top.")

    paragraphs_to_check = []
    for p in all_paragraphs[start_index:]:
        if POST_REF_PATTERN.match(p):
            if header_found or len(paragraphs_to_check) > 5:
                print(f" -> Hit trailing section boundary: '{p[:30]}...'. Stopping scan.")
                break
        paragraphs_to_check.append(p)

    pmid_pattern = re.compile(r'PMID:?\s*(\d+)', re.IGNORECASE)
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
            year_match = find_publication_year(text)
            if year_match:
                parts = extract_citation_parts(text)
                candidates = search_pubmed_by_metadata(parts)
                best_match, best_score = None, 0
                if candidates:
                    cand_meta = batch_fetch_pubmed(candidates)
                    for c_id, c_data in cand_meta.items():
                        if verify_match_precision(parts['title_snippet'], c_data['title'], parts['author'], c_data['first_author']):
                            score = fuzz.ratio(clean_for_ratio(c_data['title']), clean_for_ratio(parts['title_snippet']))
                            if score > best_score:
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
                year_match = find_publication_year(text)
                status, display_id = ("NOT_FOUND", "UNKNOWN") if year_match else ("IGNORED", display_id)
        else:
            if meta['is_retracted']: 
                status = "!! RETRACTED !!"

        if status == "IGNORED": continue 

        csv_results.append({"Original_Text": text[:60] + "...", "Final_ID": display_id, "Status": status, "Notes": "; ".join(notes)})
        
        bib_lines = [f"@article{{{item['ref_num']},"]
        if pmid: bib_lines.append(f"  pmid = {{{pmid}}},")
        if meta and meta.get('doi'): bib_lines.append(f"  doi = {{{meta['doi']}}},")
        elif item.get('found_doi'): bib_lines.append(f"  doi = {{{item['found_doi']}}},")
        
        if meta:
            clean_title = meta['title'].rstrip('.').strip()
            bib_lines.extend([f"  title = {{{clean_title}}},", f"  author = {{{meta['first_author']} et al.}},", f"  year = {{{meta['year']}}}"])
        else:
            parts = extract_citation_parts(text)
            clean_title = parts['title_snippet'].rstrip('.').strip()
            bib_lines.append(f"  title = {{{clean_title}}},")
            if parts['author']: 
                cleaned_author = clean_author_field(parts['author'])
                bib_lines.append(f"  author = {{{cleaned_author}}},")
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
    print(f"  [BibTeX] -> {bib_name}  <-- USE THIS FOR verify_bib.py")

def main():
    parser = argparse.ArgumentParser(description="Scan docx, output reports and a mapped .bib file.")
    parser.add_argument("file", help="Path to .docx", nargs='?', default=None)
    args = parser.parse_args()
    fname = args.file

    if not fname: 
        fname = input("Enter filename (.docx): ").strip().strip('"').strip("'")
    if not fname or not os.path.exists(fname): 
        sys.exit(1)

    base = re.split(r'[ _]', os.path.splitext(os.path.basename(fname))[0])[0] or "scan"
    c, t, b = analyze_document(fname)
    save_outputs(c, t, b, base)

if __name__ == "__main__":
    main()
