# RAPTOR v2.1.1 Tutorials

**Welcome to RAPTOR Tutorials!** 🦖

Step-by-step guides to master RNA-seq pipeline selection with ML-powered recommendations, ensemble analysis, **data-driven threshold optimization**, and interactive dashboards.

---

## 🆕 What's New in v2.1.1

###  Adaptive Threshold Optimizer (ATO)

v2.1.1 introduces data-driven threshold selection:
- Replace arbitrary |logFC| > 1 with optimized values
- Multiple analysis goals (discovery, balanced, validation)
- Auto-generated publication methods text
- **New Tutorial 07** dedicated to threshold optimization!

---

##  Tutorial Overview

| Tutorial | Level | Time | Focus |
|----------|-------|------|-------|
| [Tutorial 01: Getting Started](tutorial_01_getting_started.md) | Beginner | 30-45 min | ML recommendations, dashboard basics |
| [Tutorial 02: Advanced Benchmarking](tutorial_02_benchmarking.md) | Intermediate | 3-6 hours | Compare pipelines, ensemble analysis |
| [Tutorial 03: Real Data Analysis](tutorial_03_real_data.md) | Intermediate | 2-4 hours | Quality assessment, production workflows |
| [Tutorial 04: ML Recommendations](tutorial_04_ml_recommendations.md) | Advanced | 1-2 hours | Understanding ML, custom models |
| [Tutorial 05: Interactive Dashboard](tutorial_05_dashboard.md) | Intermediate | 1 hour | Web interface, visual analysis |
| [Tutorial 06: Ensemble Methods](tutorial_06_ensemble.md) | Advanced | 2-3 hours | Combining pipelines, robust results |
| [**Tutorial 07: Threshold Optimization**](tutorial_07_threshold_optimization.md) | Intermediate | 45-60 min | ** ATO - Data-driven thresholds (NEW!)** |

---

##  Learning Paths

### 🆕 **New to RAPTOR?**

1. [Tutorial 01: Getting Started](tutorial_01_getting_started.md) - Install, ML recommendations, dashboard
2. [Tutorial 07: Threshold Optimization](tutorial_07_threshold_optimization.md) ⭐ **NEW!** - Data-driven thresholds
3. [Tutorial 03: Real Data Analysis](tutorial_03_real_data.md) - Complete workflow with ATO

###  **For Publications**

1. [Tutorial 03: Real Data Analysis](tutorial_03_real_data.md) - Quality assessment
2. [Tutorial 07: Threshold Optimization](tutorial_07_threshold_optimization.md) - **Defensible thresholds**
3. [Tutorial 06: Ensemble Methods](tutorial_06_ensemble.md) - Robust results + ATO

---

##  Quick Start

### "I have DE results and need optimal thresholds" 🎯 NEW!

```python
from raptor.threshold_optimizer import optimize_thresholds
result = optimize_thresholds(de_results, goal='balanced')
print(f"Optimal: |logFC| > {result.logfc_threshold:.2f}")
print(result.methods_text)  # Copy to paper!
```

---

##  Setup

```bash
# Verify installation
raptor --version  # Should show: RAPTOR v2.1.1

# Verify ATO (NEW!)
python -c "from raptor.threshold_optimizer import optimize_thresholds; print('✅ ATO Ready!')"
```

---

##  Tutorial Version History

### v2.1.1 (Current)
-  **Added Tutorial 07: Threshold Optimization** (ATO)
-  Updated all tutorials with v2.1.1 references
-  Added ATO integration to Tutorials 01-06

---

**Happy Learning! 🦖**

*Created by Ayeh Bolouki - December 2025*
