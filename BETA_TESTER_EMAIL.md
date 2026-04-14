# Email Template — Invite Beta Testers

**Copy and customize the text below for your email/Teams/Slack message.**

---

**Subject:** Help test RAPTOR's new Data Acquisition module (20 min)

Hi [Name],

I'm looking for a few researchers to test the new **Data Acquisition module** in RAPTOR — our open-source RNA-seq analysis framework. Your feedback would be really valuable, especially since you work with [their research area].

**What it does:**
RAPTOR now lets you search GEO and SRA directly from a visual dashboard, download RNA-seq datasets with one click, pool multiple studies together with batch correction, and check data quality — all without writing code. The goal is to make cross-cohort meta-analysis accessible to any researcher.

**What I need from you:**
- 20–30 minutes to try a few scenarios (searching for studies in your field, downloading a dataset, maybe pooling two together)
- Honest feedback: what worked, what broke, what's confusing, what's missing

**How to get started:**
1. Clone the repo: `git clone https://github.com/AyehBlk/RAPTOR.git`
2. Install: `pip install -e . && pip install streamlit GEOparse biopython mygene`
3. Launch: `python -m raptor.launch_dashboard`
4. Follow the testing guide: [BETA_TESTING_GUIDE.md](https://github.com/AyehBlk/RAPTOR/blob/main/BETA_TESTING_GUIDE.md)

**How to report issues:**
- GitHub Issues (preferred): https://github.com/AyehBlk/RAPTOR/issues
- Or just reply to this email — even a one-line note is useful

**Note:** TCGA and ArrayExpress are still being built. Please focus on **GEO** and **SRA** for now.

The module is particularly designed for biomarker discovery — if you upload your own count matrix and pool it with public studies, RAPTOR handles gene ID harmonization and batch correction automatically. I'd love to hear if that workflow makes sense for your research.

Thanks!
Ayeh

---

**Shorter version (for Slack/Teams):**

Hey [Name] — I'm beta testing RAPTOR's new Data Acquisition module. It lets you search GEO/SRA, download datasets, pool studies with batch correction, all from a dashboard. Would you have 20 min to try it with your own research topic? Setup: `git clone https://github.com/AyehBlk/RAPTOR.git` then follow the testing guide in the repo. Any feedback helps — bugs, missing features, confusing UI. Thanks!
