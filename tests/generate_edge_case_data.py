"""
Generate three edge-case test datasets for RAPTOR M10 Scenario 10.

Uses the same realistic NB-count simulator as the main test data but in
three deliberately-extreme configurations:

  10a: Small sample size (5 healthy + 5 disease = 10 total, 200 genes)
       Expected: pipeline falls back to LOOCV, still delivers a panel
  10b: Many genes, few samples (200 genes × 10 samples)
       Same as 10a — the "many genes, few samples" condition IS n=10, p=200
       Expected: elastic_net picks sparse set (probably < 10)
  10c: No signal (30 healthy + 30 disease, labels shuffled randomly)
       Expected: AUC near 0.5, Youden near 0, wide bootstrap CI

10a and 10b are actually the same dataset — both test "n << p" regime and
behave the same way. We generate once.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from pathlib import Path

SEED = 42
MEAN_LIB_SIZE = 20_000_000
NB_DISPERSION = 0.15


def _nb_sample(mu, dispersion, rng):
    mu = np.maximum(mu, 1e-6)
    n = 1.0 / dispersion
    p = n / (n + mu)
    return rng.negative_binomial(n, p)


def _gene_base_expression(n_genes, rng):
    log10_tpm = rng.normal(loc=1.0, scale=0.8, size=n_genes)
    tpm = 10 ** log10_tpm
    return np.clip(tpm, 0.1, 5000)


def generate_small_cohort(seed=SEED):
    """10a/10b: 5 healthy + 5 disease, 200 genes, same signal as main data."""
    rng = np.random.default_rng(seed)
    n_healthy, n_disease = 5, 5
    n_total = n_healthy + n_disease
    n_genes = 200

    # Reuse the implanted signal from main dataset (strong + subtle)
    base_tpm = _gene_base_expression(n_genes, rng)

    conditions = ['healthy']*n_healthy + ['disease']*n_disease
    disease_mask = np.array([c == 'disease' for c in conditions])

    # Strong: BRCA1, TP53, MYC
    strong_idx = np.arange(3)
    strong_fc = np.array([+2.5, +2.2, -2.0])
    # Subtle: 8 genes
    subtle_idx = np.arange(3, 11)
    subtle_fc = rng.uniform(0.4, 1.0, 8) * rng.choice([+1, -1], 8)

    per_sample_tpm = np.tile(base_tpm, (n_total, 1))
    for gi, fc in zip(strong_idx, strong_fc):
        per_sample_tpm[disease_mask, gi] *= 2.0 ** fc
    for gi, fc in zip(subtle_idx, subtle_fc):
        per_sample_tpm[disease_mask, gi] *= 2.0 ** fc

    # Noise + library variation (scale down for small cohort)
    sample_noise = rng.normal(0, 0.40, size=(n_total, n_genes))
    per_sample_tpm *= 2.0 ** sample_noise

    lib_sizes = rng.normal(MEAN_LIB_SIZE, MEAN_LIB_SIZE * 0.4, n_total)
    lib_sizes = np.clip(lib_sizes, 0.3*MEAN_LIB_SIZE, 2.5*MEAN_LIB_SIZE)
    mu = per_sample_tpm * (lib_sizes[:, None] / 1e6)
    counts = _nb_sample(mu, NB_DISPERSION, rng)
    counts = np.clip(counts, 0, np.iinfo(np.int32).max).astype(np.int64)

    gene_names = [f"GENE_{i:03d}" for i in range(n_genes)]
    real_names = ["BRCA1","TP53","MYC","VEGFA","PTEN","RB1","CDH1","BCL2","BAX","STAT3","JAK2"]
    for i, n in enumerate(real_names):
        gene_names[i] = n

    sample_ids = [f"SMALL{i:03d}" for i in range(n_total)]
    counts_df = pd.DataFrame(counts.T, index=gene_names, columns=sample_ids)
    counts_df.index.name = "gene_id"
    meta_df = pd.DataFrame({
        "sample_id": sample_ids,
        "condition": conditions,
    })
    return counts_df, meta_df


def generate_noise_cohort(seed=SEED + 7):
    """10c: full-size cohort but labels shuffled random — no real signal."""
    rng = np.random.default_rng(seed)
    n_healthy, n_disease = 30, 30
    n_total = n_healthy + n_disease
    n_genes = 200

    base_tpm = _gene_base_expression(n_genes, rng)

    # NO condition effect applied — pure noise
    per_sample_tpm = np.tile(base_tpm, (n_total, 1))
    sample_noise = rng.normal(0, 0.40, size=(n_total, n_genes))
    per_sample_tpm *= 2.0 ** sample_noise

    lib_sizes = rng.normal(MEAN_LIB_SIZE, MEAN_LIB_SIZE * 0.4, n_total)
    lib_sizes = np.clip(lib_sizes, 0.3*MEAN_LIB_SIZE, 2.5*MEAN_LIB_SIZE)
    mu = per_sample_tpm * (lib_sizes[:, None] / 1e6)
    counts = _nb_sample(mu, NB_DISPERSION, rng)
    counts = np.clip(counts, 0, np.iinfo(np.int32).max).astype(np.int64)

    gene_names = [f"GENE_{i:03d}" for i in range(n_genes)]
    sample_ids = [f"NOISE{i:03d}" for i in range(n_total)]

    # Shuffle labels so "condition" has no real signal
    conditions = ['healthy']*n_healthy + ['disease']*n_disease
    rng.shuffle(conditions)

    counts_df = pd.DataFrame(counts.T, index=gene_names, columns=sample_ids)
    counts_df.index.name = "gene_id"
    meta_df = pd.DataFrame({
        "sample_id": sample_ids,
        "condition": conditions,
    })
    return counts_df, meta_df


# Generate all
small_counts, small_meta = generate_small_cohort()
small_counts.to_csv("test_small_counts.csv")
small_meta.to_csv("test_small_metadata.csv", index=False)

noise_counts, noise_meta = generate_noise_cohort()
noise_counts.to_csv("test_noise_counts.csv")
noise_meta.to_csv("test_noise_metadata.csv", index=False)

print("Edge-case test data generated:")
print(f"  test_small_counts.csv:    {small_counts.shape[0]} genes x {small_counts.shape[1]} samples (5h/5d)")
print(f"  test_noise_counts.csv:    {noise_counts.shape[0]} genes x {noise_counts.shape[1]} samples (shuffled, no signal)")
print(f"  Noise cohort condition balance: {noise_meta['condition'].value_counts().to_dict()}")
