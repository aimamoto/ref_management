import os
import sys
import subprocess
import argparse
import re
from pathlib import Path

def check_scripts_exist(scripts):
    """Ensure all required pipeline scripts are in the current directory."""
    missing = [s for s in scripts if not os.path.exists(s)]
    if missing:
        print("\n❌ ERROR: Missing required script(s) in the current folder:")
        for m in missing:
            print(f"  - {m}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="End-to-end AI Manuscript Reference Pipeline")
    parser.add_argument("document", type=Path, nargs='?', help="Path to the input .docx file")
    parser.add_argument("--style", choices=['1', '2', 'numbered', 'author-year'], help="1: Numbered, 2: Author-Year")
    args = parser.parse_args()

    # 1. Check Input Document
    doc_path = args.document
    if not doc_path:
        print("\n=== AI Manuscript Auto-Formatter ===")
        doc_input = input("Enter the path to your .docx file: ").strip().strip('"').strip("'")
        doc_path = Path(doc_input)
        
    if not doc_path.exists():
        print(f"\n❌ ERROR: Document '{doc_path}' not found.")
        sys.exit(1)

    # 2. Check Style Choice
    style_choice = args.style
    if not style_choice:
        print("\n=== Select Target Citation Style ===")
        print("  1. Sequential Numbered (e.g., [1-3] -> 1. Lopez-Otin C...)")
        print("  2. Author-Year (e.g., (Smith, 2023) -> Smith, I. (2023)...)")
        while True:
            choice = input("Enter 1 or 2: ").strip()
            if choice in ['1', '2']:
                style_choice = choice
                break

    # Identify active style for the final print statement
    active_style = 'author-year' if style_choice in ['2', 'author-year'] else 'numbered'

    # Required scripts
    scripts = ["scan_raw_refs_r2.py", "verify_bib_r2.py", "apply_citations_r2.py"]
    check_scripts_exist(scripts)

    # Use the current python executable to run subprocesses reliably
    python_exe = sys.executable

    print("\n" + "="*50)
    print("🚀 STARTING AUTOMATED REFERENCE PIPELINE")
    print("="*50)

    # --- STEP 1: SCAN AND EXTRACT ---
    print(f"\n>>> [1/3] Scanning {doc_path.name} for references...")
    step1 = subprocess.run([python_exe, "scan_raw_refs_r2.py", str(doc_path)])
    if step1.returncode != 0:
        print("\n❌ ERROR: Pipeline failed during Step 1 (Scanning).")
        sys.exit(1)

    # Determine the extracted output filename based on scan_raw_refs_r2.py logic
    base_name = re.split(r'[ _]', doc_path.stem)[0] or "scan"
    extracted_bib = Path(f"{base_name}_extracted.bib")
    
    if not extracted_bib.exists():
        print(f"\n❌ ERROR: Expected intermediate file '{extracted_bib}' was not generated.")
        sys.exit(1)

    # --- STEP 2: VERIFY AND ENRICH ---
    print(f"\n>>> [2/3] Fetching metadata from PubMed/Crossref for {extracted_bib.name}...")
    step2 = subprocess.run([python_exe, "verify_bib_r2.py", str(extracted_bib)])
    if step2.returncode != 0:
        print("\n❌ ERROR: Pipeline failed during Step 2 (Verification).")
        sys.exit(1)

    verified_bib = Path(f"{base_name}_extracted_verified.bib")
    if not verified_bib.exists():
         print(f"\n❌ ERROR: Expected intermediate file '{verified_bib}' was not generated.")
         sys.exit(1)

    # --- STEP 3: APPLY TO MANUSCRIPT ---
    print(f"\n>>> [3/3] Formatting Word document...")
    step3 = subprocess.run([
        python_exe, "apply_citations_r2.py", 
        str(verified_bib), 
        str(doc_path), 
        "--style", style_choice
    ])
    if step3.returncode != 0:
        print("\n❌ ERROR: Pipeline failed during Step 3 (Formatting).")
        sys.exit(1)

    # --- WRAP UP ---
    final_output = doc_path.with_name(f"{doc_path.stem}_final_{active_style}.docx")
    print("\n" + "="*50)
    print("🎉 PIPELINE COMPLETE!")
    print("="*50)
    print(f"✅ Final Document: {final_output}")
    print(f"✅ Intermediate BibTeX files ({extracted_bib}, {verified_bib}) saved for your records.")

if __name__ == "__main__":
    main()
