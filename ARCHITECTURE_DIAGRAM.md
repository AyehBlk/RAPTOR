# RAPTOR ML Recommender - System Architecture

## High-Level Overview

```
┌─────────────────────────────────────────────────────────────────────┐
│                    RAPTOR ML RECOMMENDER SYSTEM                      │
│                                                                       │
│  Intelligent Pipeline Selection Using Machine Learning               │
└─────────────────────────────────────────────────────────────────────┘

                              ┌──────────┐
                              │   USER   │
                              └────┬─────┘
                                   │
                     ┌─────────────┴─────────────┐
                     │                           │
            ┌────────▼──────────┐      ┌────────▼────────┐
            │  New RNA-seq Data │      │ Benchmark Data  │
            └────────┬──────────┘      └────────┬────────┘
                     │                           │
                     │                           │
            ┌────────▼──────────┐      ┌────────▼────────┐
            │ Data Profiling    │      │  ML Training    │
            │ (RNAseqProfiler)  │      │  Pipeline       │
            └────────┬──────────┘      └────────┬────────┘
                     │                           │
                     │                           │
            ┌────────▼──────────┐      ┌────────▼────────┐
            │ Feature Vector    │◄─────┤  Trained Model  │
            │  (30+ features)   │      │ (Random Forest) │
            └────────┬──────────┘      └─────────────────┘
                     │
                     │
            ┌────────▼──────────┐
            │   ML Predictor    │
            └────────┬──────────┘
                     │
                     │
            ┌────────▼──────────┐
            │  Recommendation   │
            │  + Confidence     │
            │  + Reasoning      │
            └───────────────────┘
```

## Detailed Component Diagram

```
┌────────────────────────────────────────────────────────────────────────┐
│                         TRAINING PHASE                                  │
└────────────────────────────────────────────────────────────────────────┘

    ┌─────────────────────┐
    │ Synthetic Generator │  or  ┌──────────────────┐
    │   (bootstrap)       │◄─────┤ Real Benchmarks  │
    └──────────┬──────────┘      └──────────────────┘
               │
               │ Generate/Collect
               ▼
    ┌──────────────────────────────────────────┐
    │      Benchmark Dataset Collection        │
    │                                          │
    │  Dataset 1: Profile + Results           │
    │  Dataset 2: Profile + Results           │
    │  Dataset 3: Profile + Results           │
    │  ...                                    │
    │  Dataset N: Profile + Results           │
    └──────────┬───────────────────────────────┘
               │
               │ Feature Extraction
               ▼
    ┌──────────────────────────────────────────┐
    │        Feature Matrix (X)                │
    │  ┌────────────────────────────────────┐ │
    │  │ sample bcv depth lib_cv zero% ... │ │
    │  │   6   0.3  high   0.15   45.2 ... │ │
    │  │   3   0.5  low    0.25   62.1 ... │ │
    │  │  12   0.2  high   0.10   38.5 ... │ │
    │  └────────────────────────────────────┘ │
    └──────────┬───────────────────────────────┘
               │
               │ Label Extraction
               ▼
    ┌──────────────────────────────────────────┐
    │        Labels (y)                        │
    │  ┌──────────────┐                       │
    │  │ best_pipeline│                       │
    │  │      3       │  (Salmon-edgeR)       │
    │  │      6       │  (NOISeq)             │
    │  │      1       │  (STAR-DESeq2)        │
    │  └──────────────┘                       │
    └──────────┬───────────────────────────────┘
               │
               │ Train-Test Split
               ▼
    ┌──────────────────────┐  ┌───────────────┐
    │   Training Set (80%) │  │ Test Set (20%)│
    └──────────┬───────────┘  └───────┬───────┘
               │                      │
               │ Train                │ Evaluate
               ▼                      │
    ┌──────────────────────┐         │
    │   Random Forest      │         │
    │   ┌──────────────┐   │         │
    │   │ Tree 1       │   │         │
    │   │ Tree 2       │   │         │
    │   │   ...        │   │         │
    │   │ Tree 200     │   │         │
    │   └──────────────┘   │         │
    └──────────┬───────────┘         │
               │                     │
               │ Validation          │
               │◄────────────────────┘
               ▼
    ┌──────────────────────┐
    │  Performance Metrics │
    │  • Accuracy: 88%     │
    │  • CV: 85% ± 2%      │
    │  • F1: 0.87          │
    └──────────┬───────────┘
               │
               │ Save
               ▼
    ┌──────────────────────┐
    │   Trained Model      │
    │   + Scaler           │
    │   + Metadata         │
    └──────────────────────┘


┌────────────────────────────────────────────────────────────────────────┐
│                        PREDICTION PHASE                                 │
└────────────────────────────────────────────────────────────────────────┘

    ┌──────────────────────┐
    │   New RNA-seq Data   │
    │  (count matrix.csv)  │
    └──────────┬───────────┘
               │
               │ Profile
               ▼
    ┌──────────────────────────────────────────┐
    │      RNAseqDataProfiler                  │
    │  ┌────────────────────────────────────┐ │
    │  │ Design: 6 samples, 20K genes      │ │
    │  │ Library: mean=20M, CV=0.15        │ │
    │  │ Counts: zero%=45, mean=1000       │ │
    │  │ BCV: 0.30                         │ │
    │  │ Depth: high (6000 reads/gene)     │ │
    │  │ Quality: 75/100                   │ │
    │  └────────────────────────────────────┘ │
    └──────────┬───────────────────────────────┘
               │
               │ Extract Features
               ▼
    ┌──────────────────────────────────────────┐
    │       Feature Extractor                  │
    │  ┌────────────────────────────────────┐ │
    │  │ n_samples: 6                      │ │
    │  │ n_genes: 20000                    │ │
    │  │ bcv: 0.30                         │ │
    │  │ depth_category: 3 (high)          │ │
    │  │ lib_size_cv: 0.15                 │ │
    │  │ zero_pct: 45.0                    │ │
    │  │ ... (24 more features)            │ │
    │  └────────────────────────────────────┘ │
    └──────────┬───────────────────────────────┘
               │
               │ Scale
               ▼
    ┌──────────────────────────────────────────┐
    │         StandardScaler                   │
    │  Normalize using training statistics     │
    └──────────┬───────────────────────────────┘
               │
               │ Predict
               ▼
    ┌──────────────────────────────────────────┐
    │        Random Forest Model               │
    │  ┌────────────────────────────────────┐ │
    │  │ Pipeline 1: 0.23                  │ │
    │  │ Pipeline 2: 0.08                  │ │
    │  │ Pipeline 3: 0.47  ◄─── Top       │ │
    │  │ Pipeline 4: 0.12                  │ │
    │  │ Pipeline 5: 0.05                  │ │
    │  │ Pipeline 6: 0.02                  │ │
    │  │ Pipeline 7: 0.02                  │ │
    │  │ Pipeline 8: 0.01                  │ │
    │  └────────────────────────────────────┘ │
    └──────────┬───────────────────────────────┘
               │
               │ Interpret
               ▼
    ┌──────────────────────────────────────────┐
    │        Recommendation Generator          │
    │  ┌────────────────────────────────────┐ │
    │  │ Primary: Pipeline 3 (47%)         │ │
    │  │   • Salmon-edgeR                  │ │
    │  │   • Fast pseudo-alignment         │ │
    │  │   • Good for high depth           │ │
    │  │                                   │ │
    │  │ Alternatives:                     │ │
    │  │   • Pipeline 1 (23%): STAR-DESeq2│ │
    │  │   • Pipeline 4 (12%): Kallisto   │ │
    │  │                                   │ │
    │  │ Key Factors:                      │ │
    │  │   • depth_category (0.156)        │ │
    │  │   • n_samples (0.132)             │ │
    │  │   • bcv (0.098)                   │ │
    │  └────────────────────────────────────┘ │
    └──────────┬───────────────────────────────┘
               │
               │ Return
               ▼
    ┌──────────────────────┐
    │   USER               │
    │ Gets recommendation  │
    │ with confidence &    │
    │ reasoning            │
    └──────────────────────┘
```

## Feature Engineering Pipeline

```
Raw Profile Data                          Feature Vector
─────────────────                        ──────────────

┌──────────────┐                         ┌──────────────┐
│ design       │──┐                      │ n_samples    │
│  n_samples=6 │  │                   ┌─▶│ n_genes      │
│  n_genes=20K │──┘                   │  │ n_conditions │
└──────────────┘                      │  │ ...          │
                                      │  └──────────────┘
┌──────────────┐                      │  
│ library_stats│──┐                   │  ┌──────────────┐
│  mean=20M    │  │                   ├─▶│ mean_lib_size│
│  cv=0.15     │──┘                   │  │ lib_size_cv  │
└──────────────┘   Feature            │  │ ...          │
                   Extractor           │  └──────────────┘
┌──────────────┐        │              │
│ count_dist   │──┐     │              │  ┌──────────────┐
│  zero%=45    │  │     │              ├─▶│ zero_pct     │
│  mean=1000   │──┘     │              │  │ mean_count   │
└──────────────┘        │              │  │ ...          │
                        │              │  └──────────────┘
┌──────────────┐        │              │
│ bio_var      │──┐     │              │  ┌──────────────┐
│  bcv=0.30    │  │     │              ├─▶│ bcv          │
│  disp=0.09   │──┘     │              │  │ dispersion   │
└──────────────┘        │              │  │ ...          │
                        ▼              │  └──────────────┘
┌──────────────┐   ┌─────────┐        │
│ sequencing   │──▶│ Extract │────────┤  ┌──────────────┐
│  depth=high  │   │  30+    │        └─▶│ total_reads  │
│  reads/g=6K  │   │Features │           │ depth_cat    │
└──────────────┘   └─────────┘           │ ...          │
                                          └──────────────┘
┌──────────────┐                          
│ complexity   │──┐                       ┌──────────────┐
│  score=75    │  │                       │ quality_score│
│  noise=0.6   │──┘                       │ noise_level  │
└──────────────┘                          │ ...          │
                                          └──────────────┘
```

## Model Decision Process

```
Input Features
      │
      │ 30+ dimensions
      ▼
┌─────────────────┐
│ Random Forest   │
│  (200 trees)    │
└────────┬────────┘
         │
         │ Each tree votes
         ▼
   ┌────────────┐    ┌────────────┐    ┌────────────┐
   │  Tree 1    │    │  Tree 2    │    │  Tree 200  │
   │            │    │            │    │            │
   │ votes: 3   │    │ votes: 3   │    │ votes: 1   │
   └────────────┘    └────────────┘    └────────────┘
         │                  │                  │
         └──────────┬───────┴──────────────────┘
                    │
                    │ Aggregate votes
                    ▼
          ┌──────────────────┐
          │ Vote Counts:     │
          │ Pipeline 1: 46   │
          │ Pipeline 2: 16   │
          │ Pipeline 3: 94 ◄──── Winner
          │ Pipeline 4: 24   │
          │ Pipeline 5: 10   │
          │ Pipeline 6:  4   │
          │ Pipeline 7:  4   │
          │ Pipeline 8:  2   │
          └────────┬─────────┘
                   │
                   │ Convert to probabilities
                   ▼
          ┌──────────────────┐
          │ Probabilities:   │
          │ Pipeline 3: 47%  │
          │ Pipeline 1: 23%  │
          │ Pipeline 4: 12%  │
          │ ... others ...   │
          └────────┬─────────┘
                   │
                   │ Add interpretation
                   ▼
          ┌──────────────────┐
          │ Recommendation   │
          └──────────────────┘
```

## Data Flow Summary

```
Training:
  Benchmarks → Features → Labels → Train → Model → Save

Prediction:
  New Data → Profile → Features → Scale → Predict → Interpret → Recommend

Integration:
  RAPTOR Profiler → ML Recommender → Enhanced CLI → User
```

## Module Dependencies

```
ml_recommender.py
  ├── numpy
  ├── pandas
  ├── scikit-learn
  │   ├── RandomForestClassifier
  │   ├── GradientBoostingClassifier
  │   └── StandardScaler
  ├── scipy
  └── joblib

synthetic_benchmarks.py
  ├── numpy
  └── pandas

example_ml_workflow.py
  ├── ml_recommender
  ├── synthetic_benchmarks
  ├── matplotlib
  └── seaborn

raptor_ml_cli.py
  ├── ml_recommender
  ├── raptor.profiler (optional)
  └── raptor.benchmark (optional)
```

## System Flow Chart

```
START
  │
  ├─► Have trained model? ──No──► Generate data ──► Train model
  │                │                                      │
  │               Yes                                     │
  │                │◄────────────────────────────────────┘
  │                │
  │                ▼
  ├─► Load count matrix
  │                │
  │                ▼
  ├─► Profile data (RNAseqDataProfiler)
  │                │
  │                ▼
  ├─► Extract features (30+ dimensions)
  │                │
  │                ▼
  ├─► Load ML model & scaler
  │                │
  │                ▼
  ├─► Predict probabilities
  │                │
  │                ▼
  ├─► Generate recommendation
  │                │
  │                ▼
  └─► Display:
      • Best pipeline (with confidence)
      • Alternatives (top 3)
      • Reasoning
      • Feature contributions
  │
END
```

This architecture diagram shows how all components work together!
