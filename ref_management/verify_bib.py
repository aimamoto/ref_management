import os
import re
import argparse
import sys
import time
import unicodedata
from pathlib import Path
from typing import Optional, Dict, List, Any, Tuple

# --- MONKEY PATCH FOR PYPARSING/BIBTEXPARSER COMPATIBILITY ---
import pyparsing
if not hasattr(pyparsing, 'DelimitedList'):
    if hasattr(pyparsing, 'delimited_list'): setattr(pyparsing, 'DelimitedList', pyparsing.delimited_list)
    elif hasattr(pyparsing, 'delimitedList'): setattr(pyparsing, 'DelimitedList', pyparsing.delimitedList)

import bibtexparser
from bibtexparser.bparser import BibTexParser
from bibtexparser.bwriter import BibTexWriter
import requests
from Bio import Entrez
from rapidfuzz import fuzz

# --- CONFIGURATION ---
Entrez.email = os.environ.get("NCBI_EMAIL", None)
Entrez.api_key = os.environ.get("NCBI_API_KEY", None)

if Entrez.api_key is None:
    print("Tip: Set NCBI_API_KEY environment variable to avoid request limits.", file=sys.stderr)

# --- HELPERS ---
def safe_entrez_call(func, *args, **kwargs):
    if not Entrez.api_key:
        time.sleep(0.35)
    else:
        time.sleep(0.1)
    return func(*args, **kwargs)

def strip_accents(s: str) -> str:
    if not s: return ""
    return ''.join(c for c in unicodedata.normalize('NFD', str(s)) if unicodedata.category(c) != 'Mn')

def clean_for_ratio(text: str) -> str:
    if not text: return ""
    text = str(text)
    text = re.sub(r'<[^>]+>', '', text)
    text = text.replace('{', '').replace('}', '')
    text = text.replace('-', ' ').replace('–', ' ').replace('/', ' ')
    text = strip_accents(text)
    return re.sub(r'[^\w\s]', '', text.lower()).strip()

def clean_author_field(author_str: str) -> str:
    if not author_str: return "Unknown"
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
        if p_clean: clean_parts.append(p_clean)
            
    if not clean_parts: return "Unknown"
    return " and ".join(clean_parts)

def get_author_surname(author_str: str) -> str:
    if not author_str: return ""
    author_str = re.sub(r'^\[?\d+\]?\s*', '', author_str.strip())
    first_block = author_str.split(',')[0].split(' and ')[0].strip()
    words = [re.sub(r'[^\w]', '', w) for w in first_block.split()]
    valid_words = [w for w in words if len(w) > 1 and w.lower() not in ['jr', 'sr', 'iii', 'ii']]
    if valid_words:
        return strip_accents(max(valid_words, key=len))
    return strip_accents(first_block[:10].strip())

def verify_match_precision(draft_title: str, fetched_title: str, draft_author: str, fetched_first_author: str) -> bool:
    title_score = fuzz.ratio(clean_for_ratio(draft_title), clean_for_ratio(fetched_title))
    if title_score < 80: return False
        
    if draft_author and fetched_first_author:
        draft_surname = get_author_surname(draft_author)
        fetched_surname = get_author_surname(fetched_first_author)
        if draft_surname and fetched_surname:
            author_score = fuzz.ratio(clean_for_ratio(draft_surname), clean_for_ratio(fetched_surname))
            if author_score < 75: return False
    return True

class BibTexEnhancer:
    def __init__(self, input_file: Path, output_file: Path):
        self.input_file = input_file
        self.output_file = output_file
        self.stats = {
            "processed": 0, "updated": 0, "failed": 0, "retracted": 0,
            "method_pmid": 0, "method_doi": 0, "method_pmc": 0,
            "method_pubmed_title": 0, "method_crossref": 0, "method_crossref_doi": 0
        }
        self.verified_titles: Dict[str, str] = {}

    def clean_text(self, text: Any) -> str:
        if not text: return ""
        text_str = str(text).replace('{', '').replace('}', '')
        return re.sub(r'\s+', ' ', text_str).strip()

    def format_pubmed_authors(self, author_list: List[str]) -> str:
        if not author_list: return ""
        formatted = []
        for name in author_list:
            parts = name.strip().split(' ')
            if len(parts) >= 2:
                surname = " ".join(parts[:-1])
                initials = parts[-1]
                formatted.append(f"{surname}, {initials}")
            else:
                formatted.append(name)
        return " and ".join(formatted)

    def search_crossref_by_doi(self, doi: str) -> Optional[Dict]:
        print(f"      ...Searching Crossref via DOI...", end='', flush=True)
        try:
            url = f"https://api.crossref.org/works/{doi}"
            resp = requests.get(url, timeout=5)
            if resp.status_code == 200:
                return resp.json().get('message', {})
        except requests.RequestException:
            pass
        return None

    def search_crossref(self, title: str, author_surname: str = "") -> Optional[List[Dict]]:
        print(f"      ...Searching Crossref...", end='', flush=True)
        try:
            url = "https://api.crossref.org/works"
            query_str = title
            if author_surname:
                query_str += f" {author_surname}"
                
            params = {
                # Upgraded to query.bibliographic which is far more reliable for exact strings than query.title alone
                "query.bibliographic": query_str.strip(), 
                "rows": 3,
                "select": "DOI,title,author,published-print,published-online,container-title,volume,issue,page"
            }
            resp = requests.get(url, params=params, timeout=5)
            if resp.status_code == 200:
                data = resp.json()
                items = data.get('message', {}).get('items', [])
                return items
        except requests.RequestException:
            pass
        return None

    def convert_crossref_to_meta(self, cr_item: Dict) -> Dict[str, str]:
        meta = {}
        meta['doi'] = cr_item.get('DOI', '')
        meta['title'] = cr_item.get('title', [''])[0]
        meta['journal'] = cr_item.get('container-title', [''])[0]

        if 'volume' in cr_item: meta['volume'] = str(cr_item['volume'])
        if 'issue' in cr_item: meta['number'] = str(cr_item['issue'])
        if 'page' in cr_item: meta['pages'] = str(cr_item['page'])

        pub = cr_item.get('published-print') or cr_item.get('published-online')
        if pub and 'date-parts' in pub: meta['year'] = str(pub['date-parts'][0][0])

        if 'author' in cr_item:
            auths = []
            for a in cr_item['author']:
                family = a.get('family', '')
                given = a.get('given', '')
                if given:
                    initial = given[0]
                    auths.append(f"{family}, {initial}")
                else:
                    auths.append(family)
            meta['author'] = " and ".join(auths)
        return meta

    def search_pubmed(self, entry: Dict[str, str]) -> Tuple[Optional[Dict], str]:
        draft_title = entry.get('title', '')

        if 'pmid' in entry and entry['pmid']:
            clean_pmid = re.sub(r'[^\d]', '', str(entry['pmid']))
            res = self.fetch_details_by_id([clean_pmid])
            if res:
                fetched_title = res.get('Title', '')
                fetched_first = res.get('AuthorList', [None])[0] if res.get('AuthorList') else ""
                if verify_match_precision(draft_title, fetched_title, entry.get('author', ''), fetched_first):
                    return (res, "pmid")

        if 'doi' in entry and entry['doi']:
            clean_doi = entry['doi'].replace('https://doi.org/', '').strip()
            res = self.fetch_by_doi(clean_doi)
            if res:
                fetched_title = res.get('Title', '')
                fetched_first = res.get('AuthorList', [None])[0] if res.get('AuthorList') else ""
                if verify_match_precision(draft_title, fetched_title, entry.get('author', ''), fetched_first):
                    return (res, "doi")
            
            cr_res = self.search_crossref_by_doi(clean_doi)
            if cr_res:
                meta = self.convert_crossref_to_meta(cr_res)
                fetched_first = meta.get('author', '').split(' and ')[0] if meta.get('author') else ""
                if verify_match_precision(draft_title, meta.get('title', ''), entry.get('author', ''), fetched_first):
                    return (meta, "crossref_doi")

        if 'url' in entry and entry['url']:
            pmc_match = re.search(r'(PMC\d+)', entry['url'])
            if pmc_match:
                res = self.fetch_by_term(pmc_match.group(1))
                if res: return (res, "pmc")

        title = self.clean_text(draft_title)
        title = re.sub(r'https?://\S+', '', title)
        
        # Fallback split to isolate title if it contains trailing journal/volume info
        title_clean_parts = re.split(r'[\.\?\!]\s+(?=[A-Z])', title)
        if title_clean_parts:
            title = title_clean_parts[0].strip()

        stopwords = {'and', 'in', 'of', 'the', 'for', 'with', 'on', 'at', 'to', 'from', 'by', 'a', 'an', 'is', 'are', 'its'}
        clean_title_str = title.replace('–', '-').replace('/', ' ')
        clean_title_str = re.sub(r'[^\w\s\-]', '', clean_title_str)
        
        # Build robust exact word-matching targeting [Title] fields to avoid literal boolean parsing clashes
        words = [w for w in clean_title_str.split() if w.lower() not in stopwords][:8]
        
        if len(words) >= 2:
            query = " AND ".join(f"{w}[Title]" for w in words)
            author = self.clean_text(entry.get('author', ''))
            surname = ""
            if author:
                surname = get_author_surname(author)
                if surname:
                    query += f" AND {surname}[Author]"
                    
            res = self.fetch_by_term(query)
            if res: return (res, "pubmed_title")

            cr_items = self.search_crossref(title, surname)
            if cr_items:
                for item in cr_items:
                    meta = self.convert_crossref_to_meta(item)
                    fetched_first = meta.get('author', '').split(' and ')[0] if meta.get('author') else ""
                    if verify_match_precision(draft_title, meta.get('title', ''), entry.get('author', ''), fetched_first):
                        return (meta, "crossref")

        return (None, "fail")

    def fetch_by_doi(self, doi: str) -> Optional[Dict]:
        try:
            handle = safe_entrez_call(Entrez.esearch, db="pubmed", term=f"{doi}[DOI]", retmax=1)
            r = Entrez.read(handle)
            if r['IdList']: return self.fetch_details_by_id(r['IdList'])
        except Exception: pass
        return None

    def fetch_by_term(self, term: str) -> Optional[Dict]:
        try:
            handle = safe_entrez_call(Entrez.esearch, db="pubmed", term=term, retmax=1)
            r = Entrez.read(handle)
            if r['IdList']: return self.fetch_details_by_id(r['IdList'])
        except Exception: pass
        return None

    def fetch_details_by_id(self, id_list: List[str]) -> Optional[Dict]:
        try:
            handle = safe_entrez_call(Entrez.esummary, db="pubmed", id=",".join(id_list), retmode="xml")
            records = Entrez.read(handle)
            return records[0] if records else None
        except Exception: return None

    def process(self):
        print(f"Reading {self.input_file}...")
        try:
            with open(self.input_file, encoding='utf-8') as bibtex_file:
                parser = BibTexParser(common_strings=True)
                bib_database = parser.parse_file(bibtex_file)
        except FileNotFoundError:
            sys.exit(1)

        total = len(bib_database.entries)
        print(f"Found {total} references. Enriching Metadata...\n")

        for index, entry in enumerate(bib_database.entries):
            entry_id = entry.get('ID', 'Unknown')
            print(f"[{index+1}/{total}] {entry_id:<25}", end='', flush=True)

            data, method = self.search_pubmed(entry)

            if data:
                self.stats[f"method_{method}"] += 1
                self.stats['updated'] += 1
                print(f" -> [OK: {method.upper()}]")

                if 'crossref' in method:
                    if data['doi']: entry['doi'] = data['doi']
                    
                    clean_title = data['title'].rstrip('.').strip()
                    clean_title = re.sub(r'<[^>]+>', '', clean_title)
                    entry['title'] = "{" + clean_title + "}"
                    
                    entry['year'] = str(data.get('year', entry.get('year', '')))
                    entry['journal'] = data.get('journal', entry.get('journal', ''))

                    if 'volume' in data: entry['volume'] = data['volume']
                    if 'number' in data: entry['number'] = data['number']
                    if 'pages' in data: entry['pages'] = data['pages']
                    if 'author' in data: entry['author'] = data['author']
                    entry['note'] = "Source: Crossref"
                else:
                    if "Retracted Publication" in data.get('PubTypeList', []):
                        print(f"      !!! WARNING: RETRACTED PUBLICATION !!!")
                        entry['note'] = f"{entry.get('note', '')} [[RETRACTED]]".strip()
                        self.stats['retracted'] += 1

                    doi_val = ""
                    if 'doi' in data.get('ArticleIds', {}):
                        doi_val = data['ArticleIds']['doi']
                    elif data.get('DOI'):
                        doi_val = data['DOI']
                        
                    if doi_val:
                        entry['doi'] = doi_val
                    entry['pmid'] = str(data.get('Id', ''))
                    
                    clean_title = data.get('Title', '').rstrip('.').strip()
                    clean_title = re.sub(r'<[^>]+>', '', clean_title)
                    entry['title'] = "{" + clean_title + "}"

                    if 'AuthorList' in data: entry['author'] = self.format_pubmed_authors(data['AuthorList'])
                    entry['journal'] = data.get('Source', entry.get('journal', ''))

                    if 'Volume' in data: entry['volume'] = str(data['Volume'])
                    if 'Issue' in data: entry['number'] = str(data['Issue'])
                    if 'Pages' in data: entry['pages'] = str(data['Pages'])

                    pub_date = data.get('PubDate', '')
                    year_match = re.search(r'\d{4}', pub_date)
                    if year_match: entry['year'] = year_match.group(0)

                clean_t = re.sub(r'[^\w]', '', entry.get('title', '').replace('{', '').replace('}', '').lower())
                if len(clean_t) > 20:
                    if clean_t in self.verified_titles:
                        print(f"      [NOTE] Exact duplicate of {self.verified_titles[clean_t]}")
                    else: self.verified_titles[clean_t] = entry_id
            else:
                self.stats['failed'] += 1
                print(" -> [NO MATCH]")

        print("\nEnrichment complete. Finalizing formatting...")
        for entry in bib_database.entries:
            author_val = entry.get('author', '').strip()
            if author_val and ',' in author_val and ' and ' not in author_val:
                entry['author'] = clean_author_field(author_val)
                
            pages_val = entry.get('pages', '').replace('{', '').replace('}', '').strip()
            if not pages_val:
                num_val = entry.get('number', '').replace('{', '').replace('}', '').strip()
                final_page = None
                
                if num_val.isdigit() and int(num_val) > 100:
                    final_page = num_val
                elif entry.get('doi'):
                    match = re.search(r'(\d+)[a-zA-Z-]*$', entry['doi'].strip())
                    if match:
                        digits = match.group(1)
                        vol = entry.get('volume', '').replace('{', '').replace('}', '').strip()
                        iss = num_val
                        if vol and digits.startswith(vol):
                            rem = digits[len(vol):]
                            if iss and rem.startswith(iss.zfill(2)):
                                rem = rem[len(iss.zfill(2)):]
                            elif iss and rem.startswith(iss):
                                rem = rem[len(iss):]
                            digits = rem if rem else digits
                        final_page = digits.lstrip('0') or digits
                
                if final_page:
                    entry['pages'] = final_page
                    entry['eid'] = final_page
                    entry['article-number'] = final_page

        print("\n" + "="*50)
        print("Completed Enrichment Verification")
        print("\n" + "="*50)
        
        writer = BibTexWriter()
        with open(self.output_file, 'w', encoding='utf-8') as f:
            f.write(writer.write(bib_database))
        print(f"Saved to: {self.output_file}")
        print(f"Verified via PubMed:   {self.stats['method_pmid'] + self.stats['method_doi'] + self.stats['method_pubmed_title']}")
        print(f"Verified via Crossref: {self.stats['method_crossref'] + self.stats['method_crossref_doi']}")
        print(f"No Match:              {self.stats['failed']}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path, nargs='?', default=None)
    args = parser.parse_args()

    if not args.output:
        args.output = args.input.with_name(f"{args.input.stem}_verified{args.input.suffix}")

    BibTexEnhancer(args.input, args.output).process()

if __name__ == "__main__":
    main()
