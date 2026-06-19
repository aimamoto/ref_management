# AI Manuscript Reference Toolkit (Revision 2)

![Python Version](https://img.shields.io/badge/python-3.x-blue) ![License](https://img.shields.io/badge/license-MIT-green)

A comprehensive Python toolkit designed to extract, verify, correct, and format references in research manuscripts. 

This **Revision 2** toolkit bridges the gap between rough drafts (which often contain raw references, incomplete metadata, or AI hallucinations) and a **finalized, submission-ready Word document**. It features an intelligent text-replacement engine that preserves your document's native fonts, italics, and formatting, while deploying "Smart Shields" to ensure scientific mathematics and amino acid codes are never accidentally converted into citation links.

## 🌟 Key Features (r2 Updates)

*   **Fully Automated Wrapper:** Run a single command to audit, enrich, and format a document end-to-end.
*   **Dual-Database Verification:** Seamlessly falls back to **Crossref** if a DOI is not found in **PubMed** (perfect for statistics or older journals).
*   **Smart Formatting Engine:** 
    *   Preserves your native Word styles, fonts, and inline italics.
    *   **Title Page Shield:** Uses NLP to skip author affiliations so they aren't confused with citations.
    *   **Math & Amino Acid Shields:** Protects $CV^2$, $R^2$, `Tyr530`, and $1 \times 10^5$ from being misread as citations.
*   **Dynamic Styles:** Interactively choose between **Sequential Numbered** (e.g., `[1–3]`) or **Author-Year** (e.g., `(Smith, 2023)`).

---

## ⚙️ Configuration (Important)

To query PubMed efficiently without hitting rate limits, you should configure your NCBI credentials.

**Option A: Environment Variables (Recommended)**
*   **Mac/Linux:**
    ```bash
    export NCBI_EMAIL="your_email@example.com"
    export NCBI_API_KEY="your_api_key"
    ```
*   **Windows (CMD/PowerShell):**
    ```cmd
    set NCBI_EMAIL=your_email@example.com
    set NCBI_API_KEY=your_api_key
    ```

**Option B: Hardcoding**
You can edit the `Entrez.email` and `Entrez.api_key` lines directly at the top of the scripts.

---

## 📦 Installation

Requires Python 3. Install the required dependencies:

```bash
pip install bibtexparser python-docx biopython rapidfuzz requests
```

---

## 🚀 Workflow 1: Fully Automated Pipeline (Recommended)

Use this wrapper script to execute the entire extraction, verification, and formatting process automatically. 

```bash
python auto_format_manuscript.py "MyDraft.docx"
```

**What it does:**
1. Prompts you to choose the citation style (Numbered or Author-Year).
2. Extracts raw references from the end of your document.
3. Downloads the missing metadata (Volume, Issue, Pages) from PubMed/Crossref.
4. Rewrites your in-text citations (e.g., collapses `[1, 2, 3]` to `[1-3]`).
5. Appends a perfectly formatted bibliography matching your document's font.

*   **Output:** `MyDraft_final_numbered.docx` (or `_author-year.docx`), plus diagnostic CSV/BibTeX files.

---

## 🔧 Workflow 2: Partial / Step-by-Step Pipeline

If you want to manually inspect or edit the references between steps, you can run the modules individually.

### Step 1: Scan & Extract
Reads the raw reference list at the bottom of your draft and maps them to PMIDs/DOIs.
```bash
python scan_raw_refs_r2.py "MyDraft.docx"
```
*   **Outputs:** 
    *   `MyDraft_extracted.bib` (Raw BibTeX generated from your document).
    *   `MyDraft_verification_report.csv` (Spreadsheet showing Retractions, Mismatches, and exact Match Scores).

### Step 2: Verify & Enrich
Takes the extracted `.bib` file, hits PubMed/Crossref, and fills in all missing Journal names, Volumes, and Authors.
```bash
python verify_bib_r2.py "MyDraft_extracted.bib"
```
*   **Output:** `MyDraft_extracted_verified.bib` *(You can open this file in a text editor to make manual tweaks if needed).*

### Step 3: Apply to Manuscript
Takes your verified references and applies them to the document.
```bash
python apply_citations_r2.py "MyDraft_extracted_verified.bib" "MyDraft.docx"
```
*   **Output:** `MyDraft_final_numbered.docx` (The final document).

---

## 📊 Extra Tools

### Generate a Text Report
If you just want a clean text file of your references (without modifying a Word document), you can use the reporter script on any verified `.bib` file:
```bash
python generate_report_r2.py "MyDraft_extracted_verified.bib"
```
*   **Output:** `MyDraft_extracted_verified_list.txt`

---

## 📂 Script Directory Overview

| Script Name | Purpose |
| :--- | :--- |
| **`auto_format_manuscript.py`** | **The Wrapper:** Runs all steps automatically. |
| **`scan_raw_refs_r2.py`** | **The Auditor:** Scans `.docx` for raw refs, outputs CSV report and a raw `.bib` mapping. |
| **`verify_bib_r2.py`** | **The Engine:** Queries PubMed/Crossref to enrich missing metadata. |
| **`apply_citations_r2.py`**| **The Formatter:** Updates inline citations, protects math/fonts, inserts Bibliography. |
| **`generate_report_r2.py`**| **The Reporter:** Converts `.bib` files into clean `.txt` lists. |

---

## 📝 Disclaimer
*While this toolkit uses fuzzy logic, NLP shields, and official APIs to verify and map data, always perform a final visual review of the generated manuscript before submitting to a journal.*
