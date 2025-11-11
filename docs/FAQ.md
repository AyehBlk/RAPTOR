# RAPTOR Frequently Asked Questions

Common questions and answers about RAPTOR.

## General Questions

**Q: What is RAPTOR?**

A: RAPTOR is a comprehensive framework for RNA-seq analysis that:
1. Implements 8 complete DE analysis pipelines
2. Provides intelligent pipeline recommendations
3. Enables systematic pipeline benchmarking
4. Makes pipeline selection evidence-based

**Q: Who should use RAPTOR?**

A: RAPTOR is for:
- Bioinformaticians analyzing RNA-seq data
- Researchers choosing analysis methods
- Method developers benchmarking approaches
- Students learning RNA-seq analysis

**Q: Is RAPTOR free?**

A: Yes! MIT license. Free and open-source.

## Installation

**Q: What are the system requirements?**

A: Minimum: 16GB RAM, 4 cores, 50GB storage  
Recommended: 32GB RAM, 8+ cores, 100GB+ storage

**Q: Do I need to install all 8 pipelines' tools?**

A: Yes, to use all features. Or install just the tools for the pipeline(s) you'll use.

**Q: Can I use Windows?**

A: Yes, via WSL2 (Windows Subsystem for Linux). Native Windows not supported for bioinformatics tools.

## Profiling & Recommendations

**Q: What input do I need?**

A: Just a count matrix (genes × samples) as CSV/TSV. Metadata optional but recommended.

**Q: Can I use FPKM or TPM values?**

A: No, must be raw integer counts. Normalized values won't work correctly.

**Q: How accurate are recommendations?**

A: Based on statistical profiling of your data + benchmark results. Typically very reliable, but always consider your specific needs.

**Q: What if I disagree with the recommendation?**

A: You can:
- Review reasoning to understand why
- Adjust scoring weights to match priorities
- Run full benchmark to compare
- Choose any pipeline you prefer

**Q: Can I use this for single-cell RNA-seq?**

A: Not yet. RAPTOR v2.0 is for bulk RNA-seq. scRNA-seq support planned for v3.0.

## Benchmarking

**Q: How long does benchmarking take?**

A: Depends on data size and pipelines:
- Quick (2 pipelines): 1-3 hours
- Full (8 pipelines): 4-24 hours

**Q: Do I need ground truth data?**

A: No, but having it allows accuracy assessment. You can benchmark without it.

**Q: Can I benchmark on simulated data?**

A: Yes! Use `raptor simulate` to generate test data.

## Pipelines

**Q: Which pipeline should I use?**

A: Use `raptor profile` to get a recommendation! But generally:
- **Most cases**: Pipeline 3 (Salmon-edgeR)
- **Highest accuracy**: Pipeline 1 (STAR-RSEM-DESeq2)
- **Large studies**: Pipeline 4 (Kallisto-Sleuth)
- **Novel transcripts**: Pipeline 2 (HISAT2-StringTie-Ballgown)

**Q: Can I add my own pipeline?**

A: Yes! See CONTRIBUTING.md for how to implement custom pipelines.

**Q: Why 8 pipelines?**

A: They represent major methodological approaches used in the field, covering alignment-based vs alignment-free, and different statistical methods.

## Technical Issues

**Q: "Command not found: raptor"**

A: Add to PATH: `export PATH=$PATH:~/.local/bin`  
Or reinstall: `pip install --force-reinstall raptor-rnaseq`

**Q: Memory errors / "Killed"**

A: Reduce resources: `raptor profile --memory 16G --threads 4`  
Or use lighter pipeline: `raptor profile --fast`

**Q: Import errors for Python packages**

A: Reinstall: `pip install --force-reinstall raptor-rnaseq[all]`

**Q: R package not found**

A: Install in R: `BiocManager::install("PackageName")`

## Data & Results

**Q: What reference genome should I use?**

A: For human: GRCh38 (hg38)  
For mouse: GRCm39 (mm39)  
Always use latest GENCODE release

**Q: How do I interpret BCV values?**

A: Biological Coefficient of Variation:
- Low (<0.2): Cell lines, controlled
- Medium (0.2-0.6): Typical studies
- High (>0.6): Clinical, variable

**Q: What's a good sequencing depth?**

A: For bulk RNA-seq:
- Minimum: 10M reads/sample
- Good: 20-30M reads/sample
- Excellent: >50M reads/sample

**Q: How many replicates do I need?**

A: Minimum 3 per group (2 in emergency)  
Recommended: 6+ per group for good power

## Citations & Publications

**Q: How do I cite RAPTOR?**

A: See CITATION.cff or:
```
Ayeh Bolouki (2025). RAPTOR: RNA-seq Analysis Pipeline 
Testing and Optimization Resource. 
GitHub: https://github.com/AyehBlk/RAPTOR
```

**Q: Can I publish results from RAPTOR?**

A: Absolutely! That's what it's for. Include:
- Which pipeline you used
- Why you chose it (recommendation + reasoning)
- Parameter settings
- RAPTOR version

## Support & Community

**Q: Where can I get help?**

A: Multiple channels:
- 📖 Documentation: docs/ folder
- 💬 GitHub Discussions
- 🐛 GitHub Issues (for bugs)
- 📧 Email: ayehbolouki1988@gmail.com

**Q: How can I contribute?**

A: See CONTRIBUTING.md! We welcome:
- Bug reports
- Feature requests
- New pipelines
- Documentation improvements
- Code contributions

**Q: Is there a mailing list?**

A: Watch the GitHub repo for updates, or follow discussions.

## Future Plans

**Q: What's coming in future versions?**

A: Roadmap:
- v2.1: Machine learning recommendations
- v3.0: Single-cell RNA-seq support
- v3.5: Long-read RNA-seq
- v4.0: Spatial transcriptomics

**Q: Can I request features?**

A: Yes! Open an issue on GitHub with your suggestion.

## Performance & Optimization

**Q: How can I speed up analysis?**

A: Tips:
- Use fast pipelines (3, 4)
- Increase threads
- Use SSD for storage
- Filter low-count genes first
- Consider subsampling for testing

**Q: How can I reduce memory usage?**

A: Options:
- Use alignment-free methods (3, 4)
- Reduce threads (counterintuitive but helps)
- Enable `memory_efficient_mode` in config
- Process samples in batches

## Miscellaneous

**Q: What does "RAPTOR" stand for?**

A: **R**NA-seq **A**nalysis **P**ipeline **T**esting and **O**ptimization **R**esource

**Q: Why a dinosaur emoji (🦖)?**

A: Raptors were smart hunters that made evidence-based decisions - just like RAPTOR helps you make informed pipeline choices!

**Q: Who maintains RAPTOR?**

A: Ayeh Bolouki (University of Namur, Belgium) with community contributions.

---

**More questions?** Check other documentation or ask on GitHub Discussions!

