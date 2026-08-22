import argparse
import re
import sys
from pathlib import Path
from typing import Dict, Any, Tuple, List

import pyparsing
if not hasattr(pyparsing, 'DelimitedList'):
    if hasattr(pyparsing, 'delimited_list'): setattr(pyparsing, 'DelimitedList', pyparsing.delimited_list)
    elif hasattr(pyparsing, 'delimitedList'): setattr(pyparsing, 'DelimitedList', pyparsing.delimitedList)

import bibtexparser
from bibtexparser.bparser import BibTexParser
from bibtexparser import customization
from rapidfuzz import fuzz

def customizations(record: Dict[str, Any]) -> Dict[str, Any]:
    record = customization.type(record)
    record = customization.author(record)
    record = customization.convert_to_unicode(record)
    return {k if k == 'ID' else k.lower(): v for k, v in record.items()}

def natural_sort_key(entry: Dict[str, Any]) -> Tuple[int, str]:
    entry_id = entry.get('ID', '')
    match = re.search(r'(\d+)', entry_id)
    return (int(match.group(1)) if match else 999999, entry_id)

def clean_latex(text: Any) -> str:
    if not text: return ""
    text_str = str(text).replace('{', '').replace('}', '')
    replacements = {r'\"a': 'ä', r'\"o': 'ö', r'\"u': 'ü', r'\%': '%', r'\&': '&', r'\_': '_'}
    for pat, rep in replacements.items(): text_str = text_str.replace(pat, rep)
    return re.sub(r'\s+', ' ', text_str.replace('\\', '')).strip()

def format_authors(entry: Dict[str, Any], style: str) -> str:
    authors = entry.get('author') or entry.get('authors')
    if not authors: return "Unknown"
    
    authors_list = authors.split(' and ') if isinstance(authors, str) else (authors if isinstance(authors, list) else [str(authors)])
    formatted_list = []
    
    for name in authors_list:
        clean_name = clean_latex(name).strip()
        if ',' in clean_name:
            parts = clean_name.split(',', 1)
            surname = parts[0].strip()
            initials = "".join([p.strip()[0].upper() for p in parts[1].split() if p.strip()])
        else:
            parts = clean_name.split()
            surname = parts[-1] if parts else "Unknown"
            initials = "".join([p.strip()[0].upper() for p in parts[:-1] if p.strip()]) if len(parts) > 1 else ""
        
        if style == 'numbered':
            formatted_list.append(f"{surname} {initials}".strip())
        else:
            formatted_list.append(f"{surname}, {'.'.join(initials)}." if initials else surname)

    if len(formatted_list) == 1: return formatted_list[0]
    elif len(formatted_list) == 2: return f"{formatted_list[0]} and {formatted_list[1]}"
    else:
        joiner = " and " if style == 'numbered' else ", and "
        return f"{', '.join(formatted_list[:-1])}{joiner}{formatted_list[-1]}"

def format_entry(entry: Dict[str, Any], style: str) -> str:
    authors = format_authors(entry, style)
    year = clean_latex(entry.get('year', '????'))
    title = clean_latex(entry.get('title', 'No Title'))
    if title and title[-1] not in '.?!': title += "."
    
    journal_str = clean_latex(entry.get('journal', ''))
    volume = clean_latex(entry.get('volume', ''))
    pages = clean_latex(entry.get('pages', '')).replace('--', '–').replace('-', '–')
    
    if volume: journal_str += f" {volume}"
    if pages: journal_str += f", {pages}" if volume else f" {pages}"
    if journal_str and journal_str[-1] not in '.?': journal_str += "."

    doi_str = f" doi: {clean_latex(entry['doi']).replace('https://doi.org/', '')}" if entry.get('doi') else ""
    
    if style == 'numbered':
        return f"{authors} ({year}) {title} {journal_str}{doi_str}".strip()
    else:
        return f"{authors} ({year}). {title} {journal_str}{doi_str}".strip()

def process_bib(input_file: Path, output_file: Path, style: str):
    print(f"Reading {input_file} (Style: {style.upper()})...")
    try:
        with open(input_file, encoding='utf-8') as f:
            parser = BibTexParser()
            parser.customization = customizations
            parser.common_strings = False
            bib = parser.parse_file(f)
    except Exception as e:
        print(f"Error parsing BibTeX: {e}"); sys.exit(1)

    # Sort based on style choice
    if style == 'numbered':
        sorted_entries = sorted(bib.entries, key=natural_sort_key)
    else:
        sorted_entries = sorted(bib.entries, key=lambda e: clean_latex(e.get('author', 'z')).lower())

    seen_pmids, seen_dois, seen_titles = {}, {}, {}

    try:
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(f"--- References List ({style.upper()} FORMAT | Source: {input_file.name}) ---\n\n")
            
            for i, entry in enumerate(sorted_entries):
                # ID tracking
                match = re.search(r'(\d+)', entry.get('ID', ''))
                ref_num = match.group(1) if match else str(i + 1)
                seq_prefix = f"{i + 1}. " if style == 'numbered' else "- "
                
                title = clean_latex(entry.get('title', '')).strip()
                if not title or (title.upper() in ['REVIEWS', 'UNKNOWN', 'REFERENCES'] and 'author' not in entry): continue

                dup_source = seen_pmids.get(entry.get('pmid')) or seen_dois.get(entry.get('doi'))
                if not dup_source and len(title) > 20:
                    norm_title = re.sub(r'[^\w]', '', title.lower())
                    for s_title, s_info in seen_titles.items():
                        if fuzz.ratio(norm_title, s_title) > 95:
                            dup_source = s_info; break
                
                if dup_source:
                    f.write(f"{seq_prefix}[DUPLICATE of {dup_source[0]}: {dup_source[1]}] {title[:50]}...\n\n")
                else:
                    if entry.get('pmid'): seen_pmids[entry['pmid']] = (ref_num, entry['ID'])
                    if entry.get('doi'): seen_dois[entry['doi']] = (ref_num, entry['ID'])
                    if len(title) > 20: seen_titles[re.sub(r'[^\w]', '', title.lower())] = (ref_num, entry['ID'])
                    
                    try:
                        f.write(f"{seq_prefix}{format_entry(entry, style)}\n\n")
                    except Exception as e:
                        f.write(f"{seq_prefix}[ERROR] {e}\n\n")

        print(f"Success! Saved to {output_file}")
    except IOError as e:
        print(f"I/O Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path, nargs='?', default=None)
    parser.add_argument("--style", choices=['1', '2', 'numbered', 'author-year'], help="1: Numbered, 2: Author-Year")
    args = parser.parse_args()
    
    if not args.input.stem.endswith("_verified"):
        candidate = args.input.with_name(f"{args.input.stem}_verified{args.input.suffix}")
        if candidate.exists(): args.input = candidate
            
    style_choice = args.style
    if not style_choice:
        print("\n=== Select Report Output Style ===")
        print("  1. Sequential Numbered (1. Lopez-Otin C...)")
        print("  2. Author-Year Alphabetical (- Smith, I. (2023)...)")
        while True:
            choice = input("Enter 1 or 2: ").strip()
            if choice in ['1', '2']:
                style_choice = choice
                break
    
    active_style = 'author-year' if style_choice in ['2', 'author-year'] else 'numbered'
    if not args.output: 
        args.output = args.input.with_name(f"{args.input.stem}_list_{active_style}.txt")
        
    process_bib(args.input, args.output, active_style)
