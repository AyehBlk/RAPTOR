#!/usr/bin/env python3
"""
RAPTOR Recommender Stress Test - 5 Scenarios
Tests rule-based and ML recommender with different RNA-seq data types.

Author: Ayeh Bolouki
"""

import sys
import json
import os
sys.path.insert(0, '.')

from raptor.recommender import PipelineRecommender, recommend_pipeline
from raptor.profiler import RNAseqDataProfiler, DataProfile

# Try ML recommender
try:
    from raptor.ml_recommender import MLPipelineRecommender
    ML_AVAILABLE = True
except ImportError:
    ML_AVAILABLE = False
    print("Note: ML recommender not available, testing rule-based only\n")

import pandas as pd
import numpy as np

os.makedirs('test_run/recommender', exist_ok=True)


def make_counts(n_genes, n_per_group, n_de=100, fold_change=2.0, 
                sparsity=0.0, lib_size_factor=1.0, seed=42):
    """Generate simulated counts with controllable properties."""
    np.random.seed(seed)
    base = np.random.gamma(shape=2, scale=200, size=n_genes)
    
    counts = {}
    for i in range(n_per_group):
        c = np.random.negative_binomial(10, 10 / (10 + base * lib_size_factor))
        counts[f'Control_{i+1}'] = c.astype(int)
    
    for i in range(n_per_group):
        modified = base.copy()
        modified[:n_de] *= fold_change
        c = np.random.negative_binomial(10, 10 / (10 + modified * lib_size_factor))
        counts[f'Treatment_{i+1}'] = c.astype(int)
    
    df = pd.DataFrame(counts, index=[f'Gene_{i+1}' for i in range(n_genes)])
    
    # Add sparsity
    if sparsity > 0:
        mask = np.random.random(df.shape) < sparsity
        df[mask] = 0
    
    return df


def make_metadata(df, confounded_batch=False):
    """Generate metadata for samples."""
    n = len(df.columns)
    half = n // 2
    conditions = ['Control'] * half + ['Treatment'] * half
    
    if confounded_batch:
        batches = ['Batch1'] * half + ['Batch2'] * half
    else:
        batches = (['Batch1', 'Batch2'] * (n // 2 + 1))[:n]
    
    return pd.DataFrame({
        'sample_id': df.columns,
        'condition': conditions,
        'batch': batches
    })


def profile_data(counts, metadata):
    """Run profiler and return DataProfile."""
    profiler = RNAseqDataProfiler(counts, metadata, group_column='condition')
    return profiler.run_full_profile()


def run_scenario(name, counts, metadata):
    """Run both recommenders on a scenario."""
    print(f"\n{'='*70}")
    print(f"  {name}")
    print(f"{'='*70}")
    print(f"  Samples: {counts.shape[1]} ({counts.shape[1]//2} per group)")
    print(f"  Genes: {counts.shape[0]:,}")
    print(f"  Zeros: {(counts == 0).sum().sum() / counts.size:.1%}")
    
    # Profile the data
    profile = profile_data(counts, metadata)
    print(f"  BCV: {profile.bcv:.3f} ({profile.bcv_category})")
    print(f"  Outliers: {profile.has_outliers}")
    print(f"  Batch effect: {profile.has_batch_effect}")
    
    # Rule-based recommendation
    recommender = PipelineRecommender(profile)
    rec = recommender.get_recommendation()
    
    print(f"\n  RULE-BASED RECOMMENDATION:")
    print(f"    Primary:     {rec.primary_pipeline} (score: {rec.primary_score:.0f})")
    print(f"    Alternative: {rec.alternative_pipeline} (score: {rec.alternative_score:.0f})")
    if rec.warnings:
        for w in rec.warnings:
            print(f"    Warning: {w}")
    
    print(f"    All scores: ", end="")
    for pipeline, score in sorted(rec.all_scores.items(), key=lambda x: -x[1]):
        print(f"{pipeline}={score:.0f}  ", end="")
    print()
    
    # ML-based recommendation
    if ML_AVAILABLE:
        try:
            ml_rec = MLPipelineRecommender()
            ml_result = ml_rec.recommend(profile)
            print(f"\n  ML RECOMMENDATION:")
            print(f"    Primary: {ml_result.primary_pipeline} (confidence: {ml_result.primary_score:.0f})")
        except Exception as e:
            print(f"\n  ML RECOMMENDATION: Error - {e}")
    
    return rec


# =============================================================================
# SCENARIO 1: Small samples (3 per group) - Should favor DESeq2
# =============================================================================
counts1 = make_counts(n_genes=2000, n_per_group=3, n_de=200, fold_change=2.0, seed=42)
meta1 = make_metadata(counts1)
rec1 = run_scenario("SCENARIO 1: Small Samples (n=3 per group)", counts1, meta1)

# =============================================================================
# SCENARIO 2: Large samples (20 per group) - Should allow limma-voom/Wilcoxon
# =============================================================================
counts2 = make_counts(n_genes=2000, n_per_group=20, n_de=200, fold_change=1.5, seed=123)
meta2 = make_metadata(counts2)
rec2 = run_scenario("SCENARIO 2: Large Samples (n=20 per group)", counts2, meta2)

# =============================================================================
# SCENARIO 3: High sparsity (60% zeros) - Should favor edgeR or Wilcoxon
# =============================================================================
counts3 = make_counts(n_genes=3000, n_per_group=6, n_de=100, sparsity=0.60, seed=456)
meta3 = make_metadata(counts3)
rec3 = run_scenario("SCENARIO 3: High Sparsity (60% zeros)", counts3, meta3)

# =============================================================================
# SCENARIO 4: Outlier present - Should favor edgeR_robust
# =============================================================================
counts4 = make_counts(n_genes=2000, n_per_group=5, n_de=150, fold_change=2.0, seed=789)
meta4 = make_metadata(counts4)
# Inject outlier
counts4['Treatment_5'] = (counts4['Treatment_5'] * 0.05).astype(int)
rec4 = run_scenario("SCENARIO 4: With Outlier Sample", counts4, meta4)

# =============================================================================
# SCENARIO 5: Batch effect confounded with condition
# =============================================================================
counts5 = make_counts(n_genes=2000, n_per_group=6, n_de=200, fold_change=2.0, seed=321)
meta5 = make_metadata(counts5, confounded_batch=True)
# Add batch effect
batch2_cols = [c for c in counts5.columns if c.startswith('Treatment')]
for col in batch2_cols:
    counts5[col] = (counts5[col] * 1.5).astype(int)
rec5 = run_scenario("SCENARIO 5: Confounded Batch Effect", counts5, meta5)


# =============================================================================
# SUMMARY
# =============================================================================
print(f"\n{'='*70}")
print(f"  SUMMARY: Recommender Results Across 5 Scenarios")
print(f"{'='*70}")
print(f"  {'Scenario':<40} {'Primary':<15} {'Score':>5}  {'Alternative':<15}")
print(f"  {'-'*80}")

scenarios = [
    ("1. Small samples (n=3)", rec1),
    ("2. Large samples (n=20)", rec2),
    ("3. High sparsity (60% zeros)", rec3),
    ("4. Outlier present", rec4),
    ("5. Confounded batch", rec5),
]

for name, rec in scenarios:
    print(f"  {name:<40} {rec.primary_pipeline:<15} {rec.primary_score:>5.0f}  {rec.alternative_pipeline:<15}")

# Expected behavior check
print(f"\n  EXPECTED BEHAVIOR CHECK:")
checks = [
    ("Small samples -> DESeq2 preferred", rec1.primary_pipeline in ['DESeq2', 'edgeR']),
    ("Large samples -> limma-voom viable", rec2.all_scores.get('limma-voom', 0) >= 80),
    ("High sparsity -> lower DESeq2 score", rec3.all_scores.get('DESeq2', 100) < rec1.all_scores.get('DESeq2', 0)),
    ("Outlier -> edgeR_robust boosted", rec4.all_scores.get('edgeR_robust', 0) > 60),
    ("Batch confounded -> warning present", any('batch' in w.lower() or 'confound' in w.lower() for w in rec5.warnings)),
]

all_pass = True
for check_name, passed in checks:
    icon = "PASS" if passed else "FAIL"
    print(f"    [{icon}] {check_name}")
    if not passed:
        all_pass = False

if all_pass:
    print(f"\n  All checks passed!")
else:
    print(f"\n  Some checks failed - recommender may need tuning")
