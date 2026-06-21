# Manuscript Reference Toolkit ARM (Another Reference Manager v1-Revision 3)

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
You can edit the `Entrez.email` and `Entrez.api_key` lines directly at the top of the `ref_management/verify_bib.py` module.

---

## 📦 Installation

Install directly from PyPI with a single command:

```bash
pip install ref-management
```

This automatically installs all required dependencies (`bibtexparser`, `python-docx`, `biopython`, `rapidfuzz`, `requests`, `citeproc-py`) and creates the following CLI commands:

| Command | Description |
| :--- | :--- |
| `arm-format` | End-to-end pipeline wrapper |
| `arm-scan` | Scan & extract references from a `.docx` |
| `arm-verify` | Enrich a `.bib` file via PubMed / Crossref |
| `arm-apply` | Apply CSL formatting to the Word document |
| `arm-add-dois` | Append missing DOIs to an intermediate draft |
| `arm-report` | Generate a plain-text reference list from a `.bib` |

---

## 🚀 Workflow 1: Fully Automated Pipeline (Recommended)

Use this wrapper command to execute the entire extraction, verification, and formatting process automatically. 

```bash
arm-format "MyDraft.docx" --csl "nature"
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
arm-scan "MyDraft.docx"
```

### Step 2: Verify & Enrich
Takes the extracted `.bib` file, hits PubMed/Crossref, and fills in all missing Journal names, Volumes, and Authors.
```bash
arm-verify "MyDraft_extracted.bib"
```

### Step 3: Apply to Manuscript via CSL Engine
Takes your verified references and applies them to the document using your target CSL style.
```bash
arm-apply "MyDraft_extracted_verified.bib" "MyDraft.docx" --csl "nature.csl"
```

---

## 🛠️ Troubleshooting

### "Dependent" CSL Style Error
If the script aborts with an error stating that your `.csl` file is a **dependent style**, it means the file you downloaded from the Zotero Style Repository is just a lightweight link to a "parent" publisher style (e.g., *The EMBO Journal* uses *EMBO Press*).

**Solution:**
1. Read the terminal error message—the script will automatically scan the XML and tell you the exact name and URL of the parent style you need.
2. Download that parent `.csl` file and place it in your `~/citation_styles/` folder.
3. Rerun the script using the parent style.

*Example:* `arm-format "MyDraft.docx" --csl embo-press`

---

## 📊 Extra Tools

### 1. Inject DOIs into an Intermediate Draft
If you want to quickly append clickable DOIs to the raw references of an intermediate draft (for co-authors to easily click/read papers) *without* fully reformatting the document or changing in-text citations:
```bash
arm-add-dois "MyDraft_extracted_verified.bib" "MyDraft.docx"
```
*   **Output:** `MyDraft_with_DOIs.docx` (Your original draft, with `https://doi.org/...` seamlessly appended to references that were missing it).

### 2. Generate a Text Report
If you just want a clean text file of your references (without modifying a Word document), you can use the reporter command on any verified `.bib` file:
```bash
arm-report "MyDraft_extracted_verified.bib"
```
*   **Output:** `MyDraft_extracted_verified_list.txt`

---

## 📂 Module Overview

| Module / Command | Purpose |
| :--- | :--- |
| **`arm-format`** | **The Wrapper:** Runs all steps automatically using the CSL Engine. |
| **`arm-scan`** | **The Auditor:** Scans `.docx` for raw refs, outputs CSV report and a raw `.bib` mapping. |
| **`arm-verify`** | **The Enrichment Engine:** Queries PubMed/Crossref to enrich missing metadata. |
| **`arm-apply`** | **The CSL Formatter:** Updates inline citations, protects math/fonts, and smartly paginates the Bibliography. |
| **`arm-add-dois`** | **The Linker:** Appends DOIs to raw reference lists for intermediate co-author drafts. |
| **`arm-report`** | **The Reporter:** Converts `.bib` files into clean `.txt` lists. |

---

## 📝 Disclaimer
*While this toolkit uses fuzzy logic, NLP shields, and official APIs to verify and map data, always perform a final visual review of the generated manuscript before submitting to a journal.*
