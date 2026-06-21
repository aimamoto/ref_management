import os
import sys
import subprocess
import argparse
import re
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Another Reference Manager (ARM)")
    parser.add_argument("document", type=Path, nargs='?', help="Path to the input .docx file")
    parser.add_argument("--csl", type=Path, help="Path to the desired .csl file")
    args = parser.parse_args()

    # 1. Check Input Document
    doc_path = args.document
    if not doc_path:
        print("\n=== ARM Manuscript Auto-Formatter (Universal CSL Edition) ===")
        doc_input = input("Enter the path to your .docx file: ").strip().strip('"').strip("'")
        doc_path = Path(doc_input)
        
    if not doc_path.exists():
        print(f"\n❌ ERROR: Document '{doc_path}' not found.")
        sys.exit(1)

    # 2. Check CSL Input
    csl_path = args.csl
    if not csl_path:
        print("\n=== Reference Formatting ===")
        csl_input = input("Enter the CSL style name or path (e.g., nature): ").strip().strip('"').strip("'")
        csl_path = Path(csl_input)

    # --- CSL Path Resolution Logic ---
    default_csl_dir = Path("~/citation_styles").expanduser()
    
    if not csl_path.exists():
        # Check if it exists in the default directory
        alt_path = default_csl_dir / csl_path.name
        if alt_path.exists():
            csl_path = alt_path
        elif not csl_path.suffix == '.csl':
            # Check if user forgot the .csl extension
            alt_path_ext = default_csl_dir / f"{csl_path.name}.csl"
            if alt_path_ext.exists():
                csl_path = alt_path_ext

    if not csl_path.exists():
        print(f"\n❌ ERROR: CSL file '{args.csl or csl_input}' not found locally or in {default_csl_dir}.")
        sys.exit(1)

    # Use the current python executable to run subprocesses reliably
    python_exe = sys.executable

    print("\n" + "="*50)
    print("🚀 STARTING AUTOMATED REFERENCE PIPELINE (CSL ENGINE)")
    print("="*50)

    # --- STEP 1: SCAN AND EXTRACT ---
    print(f"\n>>> [1/3] Scanning {doc_path.name} for references...")
    step1 = subprocess.run([python_exe, "-m", "ref_management.scan_raw_refs", str(doc_path)])
    if step1.returncode != 0:
        print("\n❌ ERROR: Pipeline failed during Step 1 (Scanning).")
        sys.exit(1)

    # Determine the extracted output filename based on scan_raw_refs logic
    base_name = re.split(r'[ _]', doc_path.stem)[0] or "scan"
    extracted_bib = Path(f"{base_name}_extracted.bib")
    
    if not extracted_bib.exists():
        print(f"\n❌ ERROR: Expected intermediate file '{extracted_bib}' was not generated.")
        sys.exit(1)

    # --- STEP 2: VERIFY AND ENRICH ---
    print(f"\n>>> [2/3] Fetching metadata from PubMed/Crossref for {extracted_bib.name}...")
    step2 = subprocess.run([python_exe, "-m", "ref_management.verify_bib", str(extracted_bib)])
    if step2.returncode != 0:
        print("\n❌ ERROR: Pipeline failed during Step 2 (Verification).")
        sys.exit(1)

    verified_bib = Path(f"{base_name}_extracted_verified.bib")
    if not verified_bib.exists():
         print(f"\n❌ ERROR: Expected intermediate file '{verified_bib}' was not generated.")
         sys.exit(1)

    # --- STEP 3: APPLY TO MANUSCRIPT USING CSL ---
    print(f"\n>>> [3/3] Formatting Word document using {csl_path.name}...")
    step3 = subprocess.run([
        python_exe, "-m", "ref_management.apply_citations",
        str(verified_bib), 
        str(doc_path), 
        "--csl", str(csl_path)
    ])
    if step3.returncode != 0:
        print("\n❌ ERROR: Pipeline failed during Step 3 (Formatting).")
        sys.exit(1)

    # --- WRAP UP ---
    final_output = doc_path.with_name(f"{doc_path.stem}_final_{csl_path.stem}.docx")
    print("\n" + "="*50)
    print("🎉 PIPELINE COMPLETE!")
    print("="*50)
    print(f"✅ Final Document: {final_output}")
    print(f"✅ Styles Applied via CSL: {csl_path.name}")
    print(f"✅ Intermediate BibTeX files ({extracted_bib}, {verified_bib}) saved for your records.")

if __name__ == "__main__":
    main()
