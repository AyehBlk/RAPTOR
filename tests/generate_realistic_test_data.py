"""
Generate realistic RNA-seq test data for RAPTOR M10 dashboard testing.

Design:
    60 samples (30 healthy, 30 disease)
    200 genes total
    - 5 STRONG DE genes   — large, easy signal (log2FC ±2.0 to ±3.0)
    - 15 SUBTLE DE genes  — weak signal (log2FC ±0.5 to ±1.5)
    - 180 NULL genes      — no signal, just biological variability

    Counts sampled from negative binomial with realistic mean-variance.
    Samples assigned to 2 batches with modest batch effects.
    Library sizes varied to mimic real sequencing yield differences.
    A few outliers injected to test QC robustness.

Validation cohort (independent 30 samples):
    Same gene identities, slightly weaker effects (replication is realistic).

Produces four files in the current directory:
    test_counts.csv
    test_metadata.csv
    test_validation_counts.csv
    test_validation_metadata.csv

Usage:
    python generate_realistic_test_data.py
"""

from __future__ import annotations
import numpy as np
import pandas as pd
from pathlib import Path

SEED = 42
rng = np.random.default_rng(SEED)

# -------- design constants ------------------------------------------------
N_HEALTHY = 30
N_DISEASE = 30
N_GENES = 200
N_STRONG = 3
N_SUBTLE = 8
N_NULL = N_GENES - N_STRONG - N_SUBTLE  # 189

# Realistic bulk-RNA-seq library size: ~20M reads per sample, ±40% variation.
MEAN_LIB_SIZE = 20_000_000
LIB_SIZE_CV = 0.40     # coefficient of variation

# Negative binomial dispersion — lower alpha = more overdispersion.
# Real bulk RNA-seq typically has alpha ~ 0.1 for most genes.
NB_DISPERSION = 0.15

# Batch effect magnitude — a subtle shift applied per batch per gene.
BATCH_EFFECT_SD = 0.40     # in log2 space — slightly stronger than before

# Per-sample biological variation in log2 space. Real bulk RNA-seq cohorts
# show sample-level log2 SD of ~0.3-0.5 across biological replicates.
SAMPLE_NOISE_SD = 0.40

# Fraction of samples to designate as outliers (library-size or expression).
OUTLIER_FRACTION = 0.05


def _nb_sample(mu: np.ndarray, dispersion: float, rng) -> np.ndarray:
    """
    Sample counts from NB(mu, dispersion) using numpy.

    Parameterization: mean = mu, variance = mu + dispersion * mu^2.
    numpy uses (n, p) where n = 1/dispersion and p = n/(n+mu).
    """
    mu = np.maximum(mu, 1e-6)
    n = 1.0 / dispersion
    p = n / (n + mu)
    return rng.negative_binomial(n, p)


def _gene_base_expression(n_genes: int, rng) -> np.ndarray:
    """
    Baseline per-gene expression (TPM-scale) drawn from a log-normal.
    Captures the fact that real genes span 5+ orders of magnitude.
    """
    log10_tpm = rng.normal(loc=1.0, scale=0.8, size=n_genes)
    tpm = 10 ** log10_tpm
    return np.clip(tpm, 0.1, 5000)


def _assign_batches(n_samples: int, rng) -> np.ndarray:
    """Split samples roughly 50/50 into two batches."""
    batch = np.array(['batch_1'] * (n_samples // 2) +
                     ['batch_2'] * (n_samples - n_samples // 2))
    rng.shuffle(batch)
    return batch


def _build_cohort(
    n_healthy: int,
    n_disease: int,
    base_tpm: np.ndarray,
    strong_idx: np.ndarray,
    subtle_idx: np.ndarray,
    strong_fc: np.ndarray,
    subtle_fc: np.ndarray,
    effect_attenuation: float,
    sample_prefix: str,
    rng,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build a (counts, metadata) pair for one cohort.

    effect_attenuation : float in (0, 1]
        Multiplier on the implanted log2FC. 1.0 = discovery cohort,
        ~0.8 = validation cohort (a mild but realistic replication dropoff).
    """
    n_samples = n_healthy + n_disease
    sample_ids = [f"{sample_prefix}{i:03d}" for i in range(n_samples)]
    conditions = ['healthy'] * n_healthy + ['disease'] * n_disease
    batches = _assign_batches(n_samples, rng)
    ages = rng.integers(25, 75, n_samples)
    sexes = rng.choice(['M', 'F'], n_samples)

    lib_sizes = rng.normal(
        loc=MEAN_LIB_SIZE,
        scale=MEAN_LIB_SIZE * LIB_SIZE_CV,
        size=n_samples,
    )
    lib_sizes = np.clip(lib_sizes, 0.3 * MEAN_LIB_SIZE, 2.5 * MEAN_LIB_SIZE)

    # A few samples get knocked down to 20-40% library size — the kind of
    # QC-borderline sample that actually appears in real GEO submissions.
    n_outliers = max(1, int(OUTLIER_FRACTION * n_samples))
    outlier_idx = rng.choice(n_samples, n_outliers, replace=False)
    lib_sizes[outlier_idx] *= rng.uniform(0.2, 0.4, n_outliers)

    # Start from baseline TPM, apply condition effects, batch effects, noise.
    per_sample_tpm = np.tile(base_tpm, (n_samples, 1))  # samples x genes

    # Condition effects: only disease samples are shifted.
    disease_mask = np.array([c == 'disease' for c in conditions])
    for gi, fc in zip(strong_idx, strong_fc):
        multiplier = 2.0 ** (fc * effect_attenuation)
        per_sample_tpm[disease_mask, gi] *= multiplier
    for gi, fc in zip(subtle_idx, subtle_fc):
        multiplier = 2.0 ** (fc * effect_attenuation)
        per_sample_tpm[disease_mask, gi] *= multiplier

    # Batch effect: additive in log2 space, gene-specific, consistent per batch.
    batch_effect = rng.normal(0, BATCH_EFFECT_SD, size=(2, N_GENES))
    batch_lookup = {'batch_1': 0, 'batch_2': 1}
    for si, b in enumerate(batches):
        per_sample_tpm[si, :] *= 2.0 ** batch_effect[batch_lookup[b], :]

    # Per-sample biological noise (sample-level jitter in log2 space).
    sample_noise = rng.normal(0, SAMPLE_NOISE_SD, size=(n_samples, N_GENES))
    per_sample_tpm *= 2.0 ** sample_noise

    # Convert TPM -> expected counts (mu) using sample-specific library sizes.
    # A TPM of 1 with a 20M library ≈ 20 expected counts for the gene.
    mu = per_sample_tpm * (lib_sizes[:, None] / 1e6)

    # Draw NB counts. Clip extreme values to int32 max for safety.
    counts = _nb_sample(mu, NB_DISPERSION, rng)
    counts = np.clip(counts, 0, np.iinfo(np.int32).max).astype(np.int64)

    # Gene names: 3 strong DE + 8 subtle DE with real symbols, rest as GENE_*.
    strong_names = ["BRCA1", "TP53", "MYC"]
    subtle_names = ["VEGFA", "PTEN", "RB1", "CDH1",
                    "BCL2", "BAX", "STAT3", "JAK2"]
    assert len(strong_names) == N_STRONG
    assert len(subtle_names) == N_SUBTLE

    gene_names = [f"GENE_{i:03d}" for i in range(N_GENES)]
    for name, gi in zip(strong_names, strong_idx):
        gene_names[gi] = name
    for name, gi in zip(subtle_names, subtle_idx):
        gene_names[gi] = name

    counts_df = pd.DataFrame(
        counts.T,  # genes x samples
        index=gene_names,
        columns=sample_ids,
    )
    counts_df.index.name = "gene_id"

    meta_df = pd.DataFrame({
        "sample_id": sample_ids,
        "condition": conditions,
        "batch": batches,
        "age": ages,
        "sex": sexes,
        "library_size": lib_sizes.astype(int),
    })

    return counts_df, meta_df, outlier_idx


# -------- gene index and effect-size assignment ---------------------------

# Reserve the first 5 gene slots for strong DE, next 15 for subtle DE.
# This keeps the gene-name assignment deterministic and simple.
strong_idx = np.arange(N_STRONG)
subtle_idx = np.arange(N_STRONG, N_STRONG + N_SUBTLE)

# Strong effects: large, mixed signs. Disease should show clear UP/DOWN.
# Pattern: 2 strong UP, 1 strong DOWN.
strong_fc = np.array([
    +2.5,   # BRCA1 up
    +2.2,   # TP53  up
    -2.0,   # MYC   down (narrative: MYC down is less usual but ok for testing)
])

# Subtle effects: log2FC drawn from a ±0.4-1.0 range, mixed signs.
# 4 UP, 4 DOWN for symmetric subtle pattern.
_subtle_magnitudes = rng.uniform(0.4, 1.0, N_SUBTLE)
_subtle_signs = np.array([+1]*4 + [-1]*4)
rng.shuffle(_subtle_signs)
subtle_fc = _subtle_magnitudes * _subtle_signs

# Base TPM — shared across cohorts so the same genes have the same baseline.
base_tpm = _gene_base_expression(N_GENES, rng)

# -------- build discovery cohort ------------------------------------------
print("Generating discovery cohort...")
disc_counts, disc_meta, disc_outliers = _build_cohort(
    n_healthy=N_HEALTHY,
    n_disease=N_DISEASE,
    base_tpm=base_tpm,
    strong_idx=strong_idx,
    subtle_idx=subtle_idx,
    strong_fc=strong_fc,
    subtle_fc=subtle_fc,
    effect_attenuation=1.0,
    sample_prefix="S",
    rng=rng,
)
disc_counts.to_csv("test_counts.csv")
disc_meta.to_csv("test_metadata.csv", index=False)

# -------- build validation cohort -----------------------------------------
print("Generating validation cohort...")
N_VAL_HEALTHY = 15
N_VAL_DISEASE = 15
val_counts, val_meta, val_outliers = _build_cohort(
    n_healthy=N_VAL_HEALTHY,
    n_disease=N_VAL_DISEASE,
    base_tpm=base_tpm,
    strong_idx=strong_idx,
    subtle_idx=subtle_idx,
    strong_fc=strong_fc,
    subtle_fc=subtle_fc,
    effect_attenuation=0.8,  # mild replication dropoff
    sample_prefix="V",
    rng=rng,
)
val_counts.to_csv("test_validation_counts.csv")
val_meta.to_csv("test_validation_metadata.csv", index=False)

# -------- print summary ---------------------------------------------------
strong_names = ["BRCA1", "TP53", "MYC"]
subtle_names = ["VEGFA", "PTEN", "RB1", "CDH1",
                "BCL2", "BAX", "STAT3", "JAK2"]

print()
print("=" * 65)
print("TEST DATA SUMMARY")
print("=" * 65)
print(f"Discovery  : {disc_counts.shape[0]} genes x {disc_counts.shape[1]} samples")
print(f"Validation : {val_counts.shape[0]} genes x {val_counts.shape[1]} samples")
print(f"Batches (discovery) : {disc_meta['batch'].value_counts().to_dict()}")
print(f"Library size (disc) : "
      f"min={disc_meta['library_size'].min():,}, "
      f"max={disc_meta['library_size'].max():,}, "
      f"median={int(disc_meta['library_size'].median()):,}")
print(f"Outlier samples (disc): {list(disc_meta['sample_id'].iloc[disc_outliers])}")
print()
print("STRONG DE genes (log2FC), expected AUC ~1.0 for these alone:")
for name, fc in zip(strong_names, strong_fc):
    direction = "UP in disease" if fc > 0 else "DOWN in disease"
    print(f"   {name:8}  log2FC={fc:+.2f}  ({direction})")
print()
print("SUBTLE DE genes (log2FC), harder to rank consistently:")
for name, fc in zip(subtle_names, subtle_fc):
    direction = "UP" if fc > 0 else "DOWN"
    print(f"   {name:8}  log2FC={fc:+.2f}  ({direction})")
print()
print(f"NULL genes : {N_NULL} (no implanted signal, biological noise only)")
print()
print("Files written:")
for f in ["test_counts.csv", "test_metadata.csv",
          "test_validation_counts.csv", "test_validation_metadata.csv"]:
    p = Path(f)
    print(f"   {f}  ({p.stat().st_size / 1024:.1f} KB)")
