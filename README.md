# Manuscript Reference Toolkit ARM (Another Reference Manager Revision 3)

![Python Version](https://img.shields.io/badge/python-3.x-blue) ![License](https://img.shields.io/badge/license-MIT-green)

A comprehensive Python toolkit designed to extract, verify, correct, and format references in research manuscripts. 

This **Revision 3** toolkit bridges the gap between rough drafts (which often contain raw references, incomplete metadata, or AI hallucinations) and a **finalized, submission-ready Word document**. It features a new **Universal CSL Formatting Engine**, an **Author-Year Bridge** (allowing you to draft with `(Author, Year)` and automatically convert to numeric formats if needed), and intelligent text-replacement algorithms that preserve your document's native fonts and formatting.

## 🌟 Key Features (r3 Updates)

*   **Universal CSL Engine:** Powered by `citeproc-py`, simply provide any Citation Style Language (`.csl`) file (e.g., from the Zotero repository) to format your manuscript exactly to specific journal requirements (e.g., Nature, Cell, APA).
*   **The Author-Year Bridge:** Draft naturally with `(Smith, 2024)` in text. The pipeline will fuzzy-match the authors to your bibliography and dynamically convert them to whatever your CSL demands (e.g., converting to `1–3` superscripts).
*   **MDPI & Online Journal Preprocessor:** Automatically algebraic-extracts article numbers from DOIs (e.g., isolating `903` from `genes15070903`) to guarantee modern online journals print with correct page numbers.
*   **Smart Bibliography Placement & Pagination:** Automatically detects trailing sections (like "Figure Legends" or "Tables") and perfectly inserts the formatted References in between the main text and trailing sections with clean page breaks.
*   **Advanced Number Collapsing:** Automatically enforces universally required typographic ranges for scientific papers (e.g., converting `1, 2, 3` into `1–3`) natively avoiding CSL engine quirks.
*   **Dual-Database Verification:** Seamlessly falls back to **Crossref** if a DOI is not found in **PubMed** (perfect for statistics or older journals).
*   **Smart Shields:** Protects $CV^2$, $R^2$, `Tyr530`, and $1 \times 10^5$ from being misread as citations.

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
You can edit the `Entrez.email` and `Entrez.api_key` lines directly at the top of the `verify_bib_r3.py` script.

---

## 📦 Installation

Requires Python 3. Install the required dependencies:

```bash
pip install bibtexparser python-docx biopython rapidfuzz requests citeproc-py
```

---

## 🚀 Workflow 1: Fully Automated Pipeline (Recommended)

Use this wrapper script to execute the entire extraction, verification, and formatting process automatically. 

```bash
python auto_format_manuscript.py "MyDraft.docx" --csl "nature"
```

> **💡 Pro-Tip (Default Directory):** You can create a folder at `~/citation_styles/` and store all your downloaded `.csl` files from Zotero there. The pipeline will automatically search this folder, meaning you can simply type `--csl cell` instead of providing the full file path.

**What it does:**
1. Loads your desired journal style via the provided `.csl` file.
2. Extracts raw references from your document.
3. Downloads the missing metadata (Volume, Issue, Pages) from PubMed/Crossref.
4. Rewrites your in-text citations natively.
5. Injects a perfectly formatted bibliography, applying proper page breaks to ensure your Tables and Figure Legends are pushed cleanly to the next page.

*   **Output:** `MyDraft_final_nature.docx`, plus diagnostic CSV/BibTeX files.

---

## 🔧 Workflow 2: Partial / Step-by-Step Pipeline

If you want to manually inspect or edit the references between steps, you can run the modules individually.

### Step 1: Scan & Extract
Reads the raw reference list at the bottom of your draft and maps them to PMIDs/DOIs.
```bash
python scan_raw_refs_r3.py "MyDraft.docx"
```

### Step 2: Verify & Enrich
Takes the extracted `.bib` file, hits PubMed/Crossref, and fills in all missing Journal names, Volumes, and Authors.
```bash
python verify_bib_r3.py "MyDraft_extract
