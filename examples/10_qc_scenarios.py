#!/usr/bin/env python3
"""
RAPTOR QC Stress Test - 4 Scenarios
Tests all QC components with realistic RNA-seq problems.

Author: Ayeh Bolouki
"""

import pandas as pd
import numpy as np
import json
import os
import sys

# Setup
sys.path.insert(0, '.')
from raptor.quality_assessment import DataQualityAssessor

def create_output(name):
    path = f'test_run/scenario_{name}'
    os.makedirs(path, exist_ok=True)
    return path


# =============================================================================
# SCENARIO 1: High-Quality Data (Gold Standard)
# 8 samples, clear condition effect, no batch effect, no outliers
# Expected: Score > 85, all components green
# =============================================================================
print("=" * 70)
print("SCENARIO 1: High-Quality RNA-seq Data")
print("=" * 70)

np.random.seed(42)
n_genes = 2000
n_per_group = 4

# Base expression with realistic distribution
base_expr = np.random.gamma(shape=2, scale=200, size=n_genes)

# Create counts with clear DE signal (200 DE genes, 2x fold change)
counts1 = {}
for i in range(n_per_group):
    counts1[f'Control_{i+1}'] = np.random.negative_binomial(
        10, 10 / (10 + base_expr)
    ).astype(int)

de_genes = np.zeros(n_genes)
de_genes[:200] = 1  # First 200 genes are DE

for i in range(n_per_group):
    modified_expr = base_expr.copy()
    modified_expr[:200] *= 2.0  # 2x fold change for DE genes
    counts1[f'Treatment_{i+1}'] = np.random.negative_binomial(
        10, 10 / (10 + modified_expr)
    ).astype(int)

df1 = pd.DataFrame(counts1, index=[f'Gene_{i+1}' for i in range(n_genes)])

# Metadata with condition only (no batch confounding)
meta1 = pd.DataFrame({
    'sample_id': df1.columns,
    'condition': ['Control']*4 + ['Treatment']*4,
    'batch': ['B1','B1','B2','B2','B1','B1','B2','B2']
})

out1 = create_output('1_high_quality')
df1.to_csv(f'{out1}/counts.csv')
meta1.to_csv(f'{out1}/metadata.csv', index=False)

assessor1 = DataQualityAssessor(df1, meta1)
result1 = assessor1.assess_quality()

print(f"  Overall Score: {result1['overall']['score']:.1f}/100")
print(f"  Status: {result1['overall']['status']}")
for name, comp in result1['components'].items():
    print(f"    {name:25s}: {comp['score']:.1f}  [{comp['status']}]")
    if comp.get('flags'):
        for f in comp['flags']:
            print(f"      -> {f}")

with open(f'{out1}/qc_results.json', 'w') as f:
    json.dump(result1, f, indent=2, default=str)
print(f"  Saved to: {out1}/")


# =============================================================================
# SCENARIO 2: Outlier Sample + Varying Library Sizes
# 1 sample has 10x lower counts, library sizes vary 5x
# Expected: Outlier detected, library quality warning
# =============================================================================
print("\n" + "=" * 70)
print("SCENARIO 2: Outlier Sample + Varying Library Sizes")
print("=" * 70)

np.random.seed(123)
df2 = df1.copy()

# Make Sample Treatment_4 a clear outlier (10x lower expression)
df2['Treatment_4'] = (df2['Treatment_4'] * 0.1).astype(int)

# Make library sizes vary dramatically
df2['Control_1'] = (df2['Control_1'] * 3.0).astype(int)  # 3x bigger
df2['Control_2'] = (df2['Control_2'] * 0.3).astype(int)  # 3x smaller

meta2 = meta1.copy()

out2 = create_output('2_outlier_libsize')
df2.to_csv(f'{out2}/counts.csv')
meta2.to_csv(f'{out2}/metadata.csv', index=False)

assessor2 = DataQualityAssessor(df2, meta2)
result2 = assessor2.assess_quality()

print(f"  Overall Score: {result2['overall']['score']:.1f}/100")
print(f"  Status: {result2['overall']['status']}")
for name, comp in result2['components'].items():
    print(f"    {name:25s}: {comp['score']:.1f}  [{comp['status']}]")
    if comp.get('flags'):
        for f in comp['flags']:
            print(f"      -> {f}")

with open(f'{out2}/qc_results.json', 'w') as f:
    json.dump(result2, f, indent=2, default=str)
print(f"  Saved to: {out2}/")


# =============================================================================
# SCENARIO 3: Strong Batch Effect (Confounded with Condition)
# Batch1 = all Controls, Batch2 = all Treatments
# Expected: Batch effect detected, confounding warning
# =============================================================================
print("\n" + "=" * 70)
print("SCENARIO 3: Confounded Batch Effect")
print("=" * 70)

np.random.seed(456)
df3 = df1.copy()

# Add strong batch effect: Batch2 samples get 1.5x expression globally
batch2_samples = ['Control_3', 'Control_4', 'Treatment_3', 'Treatment_4']
for s in batch2_samples:
    df3[s] = (df3[s] * 1.5).astype(int)

# Confounded metadata: batch perfectly correlates with condition
meta3 = pd.DataFrame({
    'sample_id': df3.columns,
    'condition': ['Control']*4 + ['Treatment']*4,
    'batch': ['Batch1']*4 + ['Batch2']*4  # Perfectly confounded!
})

out3 = create_output('3_batch_confounded')
df3.to_csv(f'{out3}/counts.csv')
meta3.to_csv(f'{out3}/metadata.csv', index=False)

assessor3 = DataQualityAssessor(df3, meta3)
result3 = assessor3.assess_quality()

print(f"  Overall Score: {result3['overall']['score']:.1f}/100")
print(f"  Status: {result3['overall']['status']}")
for name, comp in result3['components'].items():
    print(f"    {name:25s}: {comp['score']:.1f}  [{comp['status']}]")
    if comp.get('flags'):
        for f in comp['flags']:
            print(f"      -> {f}")

with open(f'{out3}/qc_results.json', 'w') as f:
    json.dump(result3, f, indent=2, default=str)
print(f"  Saved to: {out3}/")


# =============================================================================
# SCENARIO 4: Sparse Data (High Zero Inflation)
# 70% zeros, low gene detection, few samples
# Expected: Low gene detection score, warnings about sparsity
# =============================================================================
print("\n" + "=" * 70)
print("SCENARIO 4: Sparse / Low-Quality Data")
print("=" * 70)

np.random.seed(789)
n_genes4 = 3000
n_samples4 = 4

# Low expression with high zeros
base4 = np.random.gamma(shape=0.5, scale=20, size=n_genes4)
counts4 = {}
for i in range(n_samples4):
    c = np.random.negative_binomial(5, 5 / (5 + base4))
    # Add 70% zeros
    mask = np.random.random(n_genes4) < 0.70
    c[mask] = 0
    counts4[f'Sample_{i+1}'] = c.astype(int)

df4 = pd.DataFrame(counts4, index=[f'Gene_{i+1}' for i in range(n_genes4)])

meta4 = pd.DataFrame({
    'sample_id': df4.columns,
    'condition': ['A', 'A', 'B', 'B'],
    'batch': ['B1'] * 4
})

out4 = create_output('4_sparse_data')
df4.to_csv(f'{out4}/counts.csv')
meta4.to_csv(f'{out4}/metadata.csv', index=False)

assessor4 = DataQualityAssessor(df4, meta4)
result4 = assessor4.assess_quality()

print(f"  Overall Score: {result4['overall']['score']:.1f}/100")
print(f"  Status: {result4['overall']['status']}")
for name, comp in result4['components'].items():
    print(f"    {name:25s}: {comp['score']:.1f}  [{comp['status']}]")
    if comp.get('flags'):
        for f in comp['flags']:
            print(f"      -> {f}")

with open(f'{out4}/qc_results.json', 'w') as f:
    json.dump(result4, f, indent=2, default=str)
print(f"  Saved to: {out4}/")


# =============================================================================
# SUMMARY TABLE
# =============================================================================
print("\n" + "=" * 70)
print("SUMMARY: All 4 Scenarios")
print("=" * 70)
print(f"{'Scenario':<40} {'Score':>6} {'Status':<10}")
print("-" * 60)

scenarios = [
    ("1. High-Quality (Gold Standard)", result1),
    ("2. Outlier + Library Size Issues", result2),
    ("3. Confounded Batch Effect", result3),
    ("4. Sparse / Low-Quality", result4),
]

for name, r in scenarios:
    score = r['overall']['score']
    status = r['overall']['status']
    print(f"  {name:<38} {score:>5.1f}  {status}")

print("-" * 60)
print("\nAll scenario results saved to test_run/scenario_*/")
print("Each folder contains: counts.csv, metadata.csv, qc_results.json")
