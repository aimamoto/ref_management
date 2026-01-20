import bibtexparser
from bibtexparser.bparser import BibTexParser
from bibtexparser import customization
import argparse
import re
import os
from rapidfuzz import fuzz # Make sure to pip install rapidfuzz

def customizations(record):
    """Force keys to lowercase except ID."""
    record = customization.type(record)
    record = customization.author(record) 
    record = customization.convert_to_unicode(record)
    new_record = {}
    for k, v in record.items():
        if k == 'ID': new_record['ID'] = v
        else: new_record[k.lower()] = v
    return new_record

def natural_sort_key(entry):
    match = re.search(r'(\d+)', entry.get('ID', ''))
    return (int(match.group(1)) if match else 999999, entry.get('ID', ''))

def clean_latex(text):
    if not text: return ""
    text = str(text)
    while '{' in text or '}' in text:
        text = text.replace('{', '').replace('}', '')
    replacements = {r'\"a':'ä', r'\"o':'ö', r'\"u':'ü', r'\%':'%', r'\&':'&', r'\_':'_'}
    for pat, rep in replacements.items(): text = text.replace(pat, rep)
    return re.sub(r'\s+', ' ', text.replace('\\', '')).strip()

def format_authors(entry):
    authors = entry.get('author') or entry.get('authors')
    if not authors: return "Unknown"
    if isinstance(authors, str): authors = authors.split(' and ')
    
    formatted_list = []
    for name in authors:
        clean_name = clean_latex(name).strip()
        if ',' in clean_name:
            parts = clean_name.split(',', 1)
            formatted_list.append(f"{parts[0].strip()}, {parts[1].strip()[0:1].upper()}.")
        else:
            parts = clean_name.split()
            formatted_list.append(f"{parts[-1]}, {parts[0][0:1].upper()}.")

    if len(formatted_list) == 1: return formatted_list[0]
    elif len(formatted_list) == 2: return f"{formatted_list[0]} and {formatted_list[1]}"
    else: return f"{', '.join(formatted_list[:-1])}, and {formatted_list[-1]}"

def format_entry(entry):
    authors = format_authors(entry)
    year = clean_latex(entry.get('year', '????'))
    title = clean_latex(entry.get('title', 'No Title'))
    if title[-1] not in '.?': title += "."
    
    journal_str = clean_latex(entry.get('journal', ''))
    if entry.get('volume'): journal_str += f" *{clean_latex(entry['volume'])}*"
    if entry.get('pages'): journal_str += f", {clean_latex(entry['pages']).replace('--', '-')}"
    if journal_str: journal_str += "."

    ids = []
    if entry.get('pmid'): ids.append(f"PMID: {clean_latex(entry['pmid'])}")
    if entry.get('doi'): ids.append(f"doi: {clean_latex(entry['doi']).replace('https://doi.org/','')}")
    
    return f"{authors} ({year}). {title} {journal_str} {' '.join(ids)}".strip()

def process_bib(input_file, output_file):
    print(f"Reading {input_file}...")
    with open(input_file, encoding='utf-8') as f:
        parser = BibTexParser()
        parser.customization = customizations
        parser.common_strings = False
        bib = parser.parse_file(f)

    # 1. Sort Natural
    sorted_entries = sorted(bib.entries, key=natural_sort_key)
    
    # 2. Duplicate Tracking
    seen_pmids = {} # {pmid: (ref_number, ref_id)}
    seen_dois = {}
    seen_titles = {} # {clean_title: (ref_number, ref_id)}

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(f"--- References List (Source: {os.path.basename(input_file)}) ---\n\n")
        
        for i, entry in enumerate(sorted_entries):
            # Get Ref Number
            match = re.search(r'(\d+)', entry.get('ID', ''))
            ref_num = match.group(1) if match else str(i+1)
            
            # --- GARBAGE FILTER ---
            title = clean_latex(entry.get('title', '')).strip()
            if title.upper() in ['REVIEWS', 'UNKNOWN', 'REFERENCES'] and 'author' not in entry:
                continue # Skip garbage

            # --- DUPLICATE CHECK ---
            dup_source = None
            
            # Check PMID
            pmid = entry.get('pmid')
            if pmid and pmid in seen_pmids:
                dup_source = seen_pmids[pmid]
            
            # Check DOI
            doi = entry.get('doi')
            if not dup_source and doi and doi in seen_dois:
                dup_source = seen_dois[doi]

            # Check Fuzzy Title (for Preprints vs Published)
            if not dup_source:
                # Normalize title for comparison
                norm_title = re.sub(r'[^\w]', '', title.lower())
                if len(norm_title) > 20: # Only check if title is substantial
                    for s_title, s_info in seen_titles.items():
                        # High threshold (95%) to avoid false positives
                        if fuzz.ratio(norm_title, s_title) > 95:
                            dup_source = s_info
                            break
            
            # Write Entry
            if dup_source:
                f.write(f"{ref_num}. [DUPLICATE of {dup_source[0]}: {dup_source[1]}] {title[:50]}...\n\n")
            else:
                # Register this entry as seen
                if pmid: seen_pmids[pmid] = (ref_num, entry['ID'])
                if doi: seen_dois[doi] = (ref_num, entry['ID'])
                if len(title) > 20: 
                    norm_title = re.sub(r'[^\w]', '', title.lower())
                    seen_titles[norm_title] = (ref_num, entry['ID'])
                
                try:
                    line = format_entry(entry)
                    f.write(f"{ref_num}. {line}\n\n")
                except:
                    f.write(f"{ref_num}. [ERROR] {entry.get('ID')}\n\n")

    print(f"Success! Saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output", nargs='?', default=None)
    args = parser.parse_args()
    
    base, ext = os.path.splitext(args.input)
    if not base.endswith("_verified"):
        candidate = f"{base}_verified{ext}"
        if os.path.exists(candidate):
            args.input = candidate
            base, ext = os.path.splitext(candidate)
            
    if not args.output: args.output = f"{base}_list.txt"
    process_bib(args.input, args.output)
