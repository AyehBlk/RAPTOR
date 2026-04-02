# RAPTOR Beta Testing Guide — Data Acquisition Module

**Version:** 2.2.2  
**What's being tested:** GEO and SRA search, download, pooling, and quality check  
**Status:** TCGA and ArrayExpress are still under development — please focus on GEO and SRA  

---

## Thank You for Testing!

Your feedback directly shapes RAPTOR. Every bug you find and every feature you suggest makes the tool better for the entire RNA-seq community. Testing takes about 20–30 minutes.

---

## Quick Setup

### 1. Install (if you haven't already)

```bash
# Clone the repo
git clone https://github.com/AyehBlk/RAPTOR.git
cd RAPTOR

# Create environment and install
pip install -e .
pip install streamlit GEOparse biopython mygene requests pyarrow plotly
```

### 2. Launch the Dashboard

```bash
python -m raptor.launch_dashboard
# OR
python -m streamlit run raptor/dashboard/app.py
```

Navigate to the **Data Acquisition** tab (first page).

---

## Testing Scenarios

Please try as many of these as you can. For each one, note whether it **worked**, **failed**, or **could be improved**.

### Scenario 1: Search GEO for Your Research Topic

1. Go to **Search Repositories** tab
2. Select **GEO** as repository
3. Type a search query related to your research (e.g., "breast cancer RNA-seq", "Alzheimer human", "liver fibrosis mouse")
4. Click **Search**

**What to check:**
- Did results appear? How many?
- Are the results relevant to your query?
- Is the information shown (title, organism, samples, platform) useful for deciding which dataset to download?
- Is anything missing that would help you evaluate a dataset?

### Scenario 2: Download a GEO Dataset

1. From your search results, select a dataset
2. Click **Preview Info** to see study details
3. Click **Download Dataset**

**What to check:**
- Did the download complete? How long did it take?
- Does the preview show reasonable data (gene counts, sample metadata)?
- If it failed, what was the error message? Was it helpful?

### Scenario 3: Direct Accession Download

1. In the **Get Data by Accession** section at the top
2. Enter a GEO accession you already know (e.g., `GSE306761`, or any GSE from your own research)
3. Click **Preview Info**, then **Download**

**What to check:**
- Did RAPTOR recognize the accession format?
- Is the study info preview useful?
- Did the count matrix download correctly?

### Scenario 4: Search and Explore SRA

1. Switch repository to **SRA**
2. Search for a topic (e.g., "PS19 tau mouse RNA-seq")
3. Look at the results table

**What to check:**
- Did results appear?
- Is the **GEO Link** column showing linked GSE accessions?
- Select a study and click **Show Run Table** — does it load?
- Is the run information (platform, layout, read counts) useful?

### Scenario 5: Upload Your Own Data

1. Scroll down to **Upload Your Own Data**
2. Upload a count matrix CSV (genes as rows, samples as columns)
3. Optionally upload sample metadata

**What to check:**
- Did RAPTOR read your file correctly?
- Are gene counts and sample names displayed properly?
- Any issues with file format, encoding, or column detection?

### Scenario 6: Data Library Review

1. Go to the **Data Library** tab
2. Review your downloaded/uploaded datasets

**What to check:**
- Are datasets listed with correct information?
- Does the detail view (metrics, preview) make sense?
- For SRA entries: is the run table view useful?
- Does Gene ID Conversion work? (try converting between symbol/ensembl/entrez)

### Scenario 7: Pool Two or More Datasets

1. Download at least 2 GEO datasets of the **same disease/condition**
2. Go to the **Pool Datasets** tab
3. Select both datasets and click **Pool Selected Datasets**

**What to check:**
- Did pooling complete?
- Does the gene overlap preview make sense?
- Did batch correction options work?
- Any errors during pooling?

### Scenario 8: Quality Check

1. After pooling, go to the **Quality Check** tab
2. Review the QC plots (library sizes, PCA, correlations, distributions)

**What to check:**
- Do the plots render correctly?
- Is the quality verdict helpful?
- Can you see batch effects in the PCA?
- Any suggestions for additional QC metrics?

### Scenario 9: Export

1. Go to the **Export** tab
2. Try downloading counts as CSV
3. Try storing for downstream analysis

**What to check:**
- Do downloads work?
- Is the file format correct when you open it in Excel/R?

---

## What to Report

For **bugs** (something broke):
- What you were trying to do
- What happened instead
- The exact error message (screenshot or copy-paste)
- Your dataset accession (e.g., GSE306761) so we can reproduce it
- Your OS (Windows/Mac/Linux)

For **feature requests** (something that would be useful):
- What you were trying to accomplish
- Why the current interface doesn't support it well
- How you'd like it to work

For **general impressions**:
- Was anything confusing?
- What did you like?
- What would make your meta-analysis workflow easier?

---

## How to Report

**Option A — GitHub Issues (preferred):**  
Go to [github.com/AyehBlk/RAPTOR/issues](https://github.com/AyehBlk/RAPTOR/issues) and click **New Issue**. Choose "Bug Report" or "Feature Request" — the template will guide you.

**Option B — Email:**  
Send your feedback to ayehbolouki1988@gmail.com with subject line: **RAPTOR Beta Feedback**

**Option C — Quick note:**  
Even a one-line message like "GSE123456 download failed with timeout error" is valuable. Don't feel you need to write a full report.

---

## Known Limitations (Don't Report These)

- **TCGA and ArrayExpress** connectors are still under development — they may not work reliably yet
- **Gene ID conversion** requires the `mygene` package (`pip install mygene`)
- **GEO download** can be slow for very large datasets (>100 samples) — this is normal
- **SRA studies** show run metadata, not count matrices — this is by design (the linked GSE has the counts)
- The `raptor dashboard` CLI command may not work — use `python -m raptor.launch_dashboard` instead

---

## Your Feedback Matters

Every piece of feedback helps. The goal is to make RAPTOR's Data Acquisition module flexible enough to handle any RNA-seq meta-analysis workflow — from simple single-study downloads to complex multi-cohort biomarker discovery. Your real-world use cases are the best test.

Thank you!

— Ayeh
