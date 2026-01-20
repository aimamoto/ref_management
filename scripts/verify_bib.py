import bibtexparser
from bibtexparser.bwriter import BibTexWriter
import os
import re
import argparse
import requests
import io
from Bio import Entrez
from pypdf import PdfReader
from rapidfuzz import fuzz

# --- CONFIGURATION ---
Entrez.email = os.environ.get("NCBI_EMAIL", "aimamoto@uchicago.edu")
Entrez.api_key = os.environ.get("NCBI_API_KEY", None)

if Entrez.api_key is None:
    print("Tip: Set NCBI_API_KEY environment variable to avoid request limits.")

# --- HELPERS ---
class BibTexEnhancer:
    def __init__(self, input_file, output_file):
        self.input_file = input_file
        self.output_file = output_file
        self.stats = {
            "processed": 0, "updated": 0, "failed": 0, "retracted": 0,
            "method_pmid": 0, "method_doi": 0, "method_pmc": 0,
            "method_pdf": 0, "method_pubmed_title": 0, "method_crossref": 0
        }
        self.verified_titles = {} 

    def clean_text(self, text):
        if not text: return ""
        text = str(text)
        text = text.replace('{', '').replace('}', '')
        return re.sub(r'\s+', ' ', text).strip()

    def format_pubmed_authors(self, author_list):
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

    def search_crossref(self, title):
        """Searches Crossref API for non-PubMed papers."""
        print(f"      ...Searching Crossref...", end='', flush=True)
        try:
            url = "https://api.crossref.org/works"
            params = {
                "query.title": title,
                "rows": 1,
                "select": "DOI,title,author,published-print,published-online,container-title,volume,issue,page"
            }
            resp = requests.get(url, params=params, timeout=5)
            if resp.status_code == 200:
                data = resp.json()
                items = data.get('message', {}).get('items', [])
                if items:
                    item = items[0]
                    found_title = item.get('title', [''])[0]
                    if fuzz.ratio(title.lower(), found_title.lower()) > 80:
                        return item
        except: pass
        return None

    def convert_crossref_to_meta(self, cr_item):
        """Converts Crossref JSON to our internal dict format."""
        meta = {}
        meta['doi'] = cr_item.get('DOI', '')
        meta['title'] = cr_item.get('title', [''])[0]
        meta['journal'] = cr_item.get('container-title', [''])[0]
        
        # Bibliographic Details
        if 'volume' in cr_item: meta['volume'] = str(cr_item['volume'])
        if 'issue' in cr_item: meta['number'] = str(cr_item['issue'])
        if 'page' in cr_item: meta['pages'] = str(cr_item['page'])

        # Year
        pub = cr_item.get('published-print') or cr_item.get('published-online')
        if pub and 'date-parts' in pub:
            meta['year'] = str(pub['date-parts'][0][0])
        
        # Authors
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

    def search_pubmed(self, entry):
        # 1. Existing PMID
        if 'pmid' in entry and entry['pmid']:
            clean_pmid = re.sub(r'[^\d]', '', str(entry['pmid']))
            res = self.fetch_details_by_id([clean_pmid])
            if res: return (res, "pmid")

        # 2. Existing DOI
        if 'doi' in entry and entry['doi']:
            clean_doi = entry['doi'].replace('https://doi.org/', '').strip()
            res = self.fetch_by_doi(clean_doi)
            if res: return (res, "doi")

        # 3. PMC URL
        if 'url' in entry and entry['url']:
            pmc_match = re.search(r'(PMC\d+)', entry['url'])
            if pmc_match:
                res = self.fetch_by_term(pmc_match.group(1))
                if res: return (res, "pmc")

        # 4. Search by Title (PubMed)
        title = self.clean_text(entry.get('title', ''))
        search_title = re.sub(r'[^\w\s]', '', title)
        
        if len(search_title) > 10:
            query = f"{search_title}[Title]"
            author = self.clean_text(entry.get('author', ''))
            if author:
                first_author = author.split(',')[0].split(' ')[-1]
                query += f" AND {first_author}[Author]"
            
            res = self.fetch_by_term(query)
            if res: return (res, "pubmed_title")

        # 5. Search by Title (Crossref)
        if len(search_title) > 10:
            cr_res = self.search_crossref(search_title)
            if cr_res:
                return (self.convert_crossref_to_meta(cr_res), "crossref")

        return (None, "fail")

    def fetch_by_doi(self, doi):
        try:
            handle = Entrez.esearch(db="pubmed", term=f"{doi}[DOI]", retmax=1)
            r = Entrez.read(handle)
            if r['IdList']: return self.fetch_details_by_id(r['IdList'])
        except: pass
        return None

    def fetch_by_term(self, term):
        try:
            handle = Entrez.esearch(db="pubmed", term=term, retmax=1)
            r = Entrez.read(handle)
            if r['IdList']: return self.fetch_details_by_id(r['IdList'])
        except: pass
        return None

    def fetch_details_by_id(self, id_list):
        try:
            handle = Entrez.esummary(db="pubmed", id=",".join(id_list), retmode="xml")
            records = Entrez.read(handle)
            return records[0] if records else None
        except: return None

    def process(self):
        print(f"Reading {self.input_file}...")
        with open(self.input_file, encoding='utf-8') as bibtex_file:
            parser = bibtexparser.bparser.BibTexParser(common_strings=True)
            bib_database = parser.parse_file(bibtex_file)

        total = len(bib_database.entries)
        print(f"Found {total} references. Enriching Metadata...\n")

        for index, entry in enumerate(bib_database.entries):
            entry_id = entry.get('ID', 'Unknown')
            print(f"[{index+1}/{total}] {entry_id:<25}", end='', flush=True)
            
            data, method = self.search_pubmed(entry)

            if data:
                # Update Stats
                self.stats[f"method_{method}"] += 1
                self.stats['updated'] += 1
                print(f" -> [OK: {method.upper()}]")

                if method == 'crossref':
                    # --- Crossref Mapping ---
                    if data['doi']: entry['doi'] = data['doi']
                    entry['title'] = "{" + data['title'] + "}"
                    entry['year'] = str(data.get('year', entry.get('year', '')))
                    entry['journal'] = data.get('journal', entry.get('journal', ''))
                    
                    if 'volume' in data: entry['volume'] = data['volume']
                    if 'number' in data: entry['number'] = data['number']
                    if 'pages' in data: entry['pages'] = data['pages']
                    if 'author' in data: entry['author'] = data['author']
                    entry['note'] = "Source: Crossref"
                else:
                    # --- PubMed Mapping ---
                    if "Retracted Publication" in data.get('PubTypeList', []):
                        print(f"      !!! WARNING: RETRACTED PUBLICATION !!!")
                        entry['note'] = f"{entry.get('note', '')} [[RETRACTED]]".strip()
                        self.stats['retracted'] += 1

                    if 'doi' in data.get('ArticleIds', {}):
                        entry['doi'] = data['ArticleIds']['doi']
                    entry['pmid'] = str(data.get('Id', ''))
                    entry['title'] = "{" + data.get('Title', '') + "}"
                    
                    if 'AuthorList' in data:
                        entry['author'] = self.format_pubmed_authors(data['AuthorList'])
                    entry['journal'] = data.get('Source', entry.get('journal', ''))
                    
                    # --- NEW: Extract Vol/Issue/Pages from PubMed XML ---
                    if 'Volume' in data: entry['volume'] = str(data['Volume'])
                    if 'Issue' in data: entry['number'] = str(data['Issue']) # BibTeX uses 'number' for issue
                    if 'Pages' in data: entry['pages'] = str(data['Pages'])
                    
                    pub_date = data.get('PubDate', '')
                    year_match = re.search(r'\d{4}', pub_date)
                    if year_match: entry['year'] = year_match.group(0)

                # Deduplication check
                clean_t = re.sub(r'[^\w]', '', entry.get('title', '').replace('{','').replace('}','').lower())
                if len(clean_t) > 20:
                    if clean_t in self.verified_titles:
                        orig_id = self.verified_titles[clean_t]
                        print(f"      [NOTE] Exact duplicate of {orig_id}")
                    else:
                        self.verified_titles[clean_t] = entry_id

            else:
                self.stats['failed'] += 1
                print(" -> [NO MATCH]")

        print("\n" + "="*50)
        writer = BibTexWriter()
        with open(self.output_file, 'w', encoding='utf-8') as bibfile:
            bibfile.write(writer.write(bib_database))

        print(f"Saved to: {self.output_file}")
        print(f"Verified via PubMed:   {self.stats['method_pmid'] + self.stats['method_doi'] + self.stats['method_pubmed_title']}")
        print(f"Verified via Crossref: {self.stats['method_crossref']}")
        print(f"No Match:              {self.stats['failed']}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output", nargs='?', default=None)
    args = parser.parse_args()
    
    if not args.output:
        base, ext = os.path.splitext(args.input)
        args.output = f"{base}_verified{ext}"

    BibTexEnhancer(args.input, args.output).process()
