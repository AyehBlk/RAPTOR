"""
RAPTOR Module 9: Ensemble Analysis
===================================

This example demonstrates how to combine differential expression results
from multiple methods using ensemble techniques to create high-confidence
consensus gene lists.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0

What this example covers:
- Fisher's method (p-value combination)
- Brown's method (correlation-aware)
- Robust Rank Aggregation (RRA)
- Simple Voting (count-based)
- Weighted Ensemble (performance-based)
- Method comparison and selection
- Creating validation priority lists
"""

import pandas as pd
import pickle
import numpy as np
from pathlib import Path

from raptor.ensemble import (
    ensemble_fisher,
    ensemble_brown,
    ensemble_rra,
    ensemble_voting,
    ensemble_weighted
)

print("=" * 70)
print("RAPTOR Module 9: Ensemble Analysis Example")
print("=" * 70)
print()

# =============================================================================
# PART 1: Load DE Results
# =============================================================================

print("PART 1: Loading DE Results")
print("-" * 70)
print()

# Load DE results from Module 7
print("Loading DE results...")

with open('output/module7_imported/deseq2_result.pkl', 'rb') as f:
    deseq2_result = pickle.load(f)

with open('output/module7_imported/edger_result.pkl', 'rb') as f:
    edger_result = pickle.load(f)

with open('output/module7_imported/limma_result.pkl', 'rb') as f:
    limma_result = pickle.load(f)

# Organize into dictionary
de_results = {
    'DESeq2': deseq2_result,
    'edgeR': edger_result,
    'limma': limma_result
}

print(f"✓ Loaded results from {len(de_results)} methods")
print()

# Quick summary of individual methods
for method, result in de_results.items():
    n_sig = len(result.data[result.data['padj'] < 0.05])
    print(f"  {method}: {n_sig} significant genes (FDR < 0.05)")
print()

# =============================================================================
# PART 2: Fisher's Method (Maximum Sensitivity)
# =============================================================================

print("PART 2: Fisher's Method - P-value Combination")
print("-" * 70)
print()

print("Fisher's method combines p-values from multiple tests.")
print("Best for: Maximum sensitivity, exploratory analysis")
print()

fisher_result = ensemble_fisher(
    de_results=de_results,
    use_padj=False,              # IMPORTANT: Use raw p-values
    significance_threshold=0.05,
    check_direction=True,
    direction_threshold=1.0,     # All methods must agree
    output_dir='output/module9_ensemble/fisher/'
)

print(f"✓ Fisher's method found {fisher_result.n_consensus_genes} genes")
print()

# =============================================================================
# PART 3: Brown's Method (Correlation-Aware)
# =============================================================================

print("PART 3: Brown's Method - Correlation-Aware Combination")
print("-" * 70)
print()

print("Brown's method accounts for correlation between methods.")
print("Best for: When methods use the same data/normalization")
print()

brown_result = ensemble_brown(
    de_results=de_results,
    use_padj=False,
    significance_threshold=0.05,
    check_direction=True,
    direction_threshold=1.0,
    output_dir='output/module9_ensemble/brown/'
)

print(f"✓ Brown's method found {brown_result.n_consensus_genes} genes")
print()

# =============================================================================
# PART 4: Robust Rank Aggregation (RRA)
# =============================================================================

print("PART 4: Robust Rank Aggregation (RRA)")
print("-" * 70)
print()

print("RRA combines ranked gene lists using order statistics.")
print("Best for: Robust ranking, handling outliers")
print()

rra_result = ensemble_rra(
    de_results=de_results,
    rank_by='pvalue',
    use_padj=False,
    significance_threshold=0.05,
    check_direction=True,
    output_dir='output/module9_ensemble/rra/'
)

print(f"✓ RRA found {rra_result.n_consensus_genes} genes")
print()

# =============================================================================
# PART 5: Simple Voting (High Confidence)
# =============================================================================

print("PART 5: Simple Voting - Count-Based Consensus")
print("-" * 70)
print()

# Strategy 1: Unanimous (all 3 methods)
print("Strategy 1: Unanimous voting (all 3 methods must agree)")

voting_unanimous = ensemble_voting(
    de_results=de_results,
    min_methods=3,
    filters={'padj': 0.05, 'lfc': 0.0},
    check_direction=True,
    direction_threshold=1.0,
    output_dir='output/module9_ensemble/voting_unanimous/'
)

print(f"  ✓ Unanimous: {voting_unanimous.n_consensus_genes} genes")

# Strategy 2: Majority (≥2 methods)
print("Strategy 2: Majority voting (≥2 methods)")

voting_majority = ensemble_voting(
    de_results=de_results,
    min_methods=2,
    filters={'padj': 0.05, 'lfc': 0.0},
    check_direction=True,
    direction_threshold=0.67,
    output_dir='output/module9_ensemble/voting_majority/'
)

print(f"  ✓ Majority: {voting_majority.n_consensus_genes} genes")
print()

# =============================================================================
# PART 6: Weighted Ensemble (Performance-Based)
# =============================================================================

print("PART 6: Weighted Ensemble - Performance-Based")
print("-" * 70)
print()

# Try to load weights from Module 8
try:
    weights_df = pd.read_csv('output/module8_optimization/ensemble_weights.csv')
    weights = dict(zip(weights_df['Method'], weights_df['Weight']))
    print("Using weights from Module 8 optimization:")
    for method, weight in weights.items():
        print(f"  {method}: {weight:.4f}")
    print()
except FileNotFoundError:
    print("No Module 8 weights found, using equal weights")
    weights = None
    print()

weighted_result = ensemble_weighted(
    de_results=de_results,
    weights=weights,
    min_score=1.5,
    filters={'padj': 0.05},
    check_direction=True,
    output_dir='output/module9_ensemble/weighted/'
)

print(f"✓ Weighted ensemble found {weighted_result.n_consensus_genes} genes")
print()

# =============================================================================
# PART 7: Compare All Methods
# =============================================================================

print("PART 7: Method Comparison")
print("-" * 70)
print()

# Organize results
ensemble_results = {
    'Fisher': fisher_result,
    'Brown': brown_result,
    'RRA': rra_result,
    'Voting (≥2)': voting_majority,
    'Voting (≥3)': voting_unanimous,
    'Weighted': weighted_result
}

# Create comparison table
comparison_data = []
for method_name, result in ensemble_results.items():
    stats = result.method_statistics.get('overall', {})
    comparison_data.append({
        'Method': method_name,
        'Genes': result.n_consensus_genes,
        'Upregulated': stats.get('n_upregulated', 0),
        'Downregulated': stats.get('n_downregulated', 0)
    })

comparison_df = pd.DataFrame(comparison_data)

print("Performance Comparison:")
print()
print(comparison_df.to_string(index=False))
print()

# =============================================================================
# PART 8: Overlap Analysis
# =============================================================================

print("PART 8: Overlap Analysis")
print("-" * 70)
print()

# Get gene sets
gene_sets = {
    'Fisher': set(fisher_result.consensus_genes['gene_id']),
    'Brown': set(brown_result.consensus_genes['gene_id']),
    'RRA': set(rra_result.consensus_genes['gene_id']),
    'Voting_2': set(voting_majority.consensus_genes['gene_id']),
    'Voting_3': set(voting_unanimous.consensus_genes['gene_id']),
    'Weighted': set(weighted_result.consensus_genes['gene_id'])
}

# Core consensus (all methods agree)
core_consensus = set.intersection(*gene_sets.values())
print(f"Core consensus (all 6 strategies): {len(core_consensus)} genes")
print()

# Pairwise overlaps
print("Pairwise overlaps (Jaccard index):")
print()

methods = list(gene_sets.keys())
for i, method1 in enumerate(methods):
    for method2 in methods[i+1:]:
        overlap = len(gene_sets[method1] & gene_sets[method2])
        union = len(gene_sets[method1] | gene_sets[method2])
        jaccard = overlap / union if union > 0 else 0
        print(f"  {method1:12s} ∩ {method2:12s}: {overlap:4d} genes (J={jaccard:.3f})")

print()

# =============================================================================
# PART 9: Create Tiered Gene List
# =============================================================================

print("PART 9: Creating Tiered Gene List")
print("-" * 70)
print()

print("Creating 3-tier consensus gene list...")
print()

# Tier 1: Core consensus (all 6 strategies)
tier1_genes = list(core_consensus)
tier1_df = pd.DataFrame({'gene_id': tier1_genes})
tier1_df['tier'] = 1
tier1_df['confidence'] = 'Very High'
tier1_df['n_methods'] = 6
print(f"  Tier 1 (Very High): {len(tier1_genes)} genes (all 6 agree)")

# Tier 2: Detected by ≥4 strategies
tier2_genes = []
for gene in set.union(*gene_sets.values()):
    count = sum(1 for gs in gene_sets.values() if gene in gs)
    if 4 <= count < 6:
        tier2_genes.append(gene)

tier2_df = pd.DataFrame({'gene_id': tier2_genes})
tier2_df['tier'] = 2
tier2_df['confidence'] = 'High'
tier2_df['n_methods'] = tier2_df['gene_id'].apply(
    lambda g: sum(1 for gs in gene_sets.values() if g in gs)
)
print(f"  Tier 2 (High): {len(tier2_genes)} genes (4-5 agree)")

# Tier 3: Detected by 2-3 strategies
tier3_genes = []
for gene in set.union(*gene_sets.values()):
    count = sum(1 for gs in gene_sets.values() if gene in gs)
    if 2 <= count < 4:
        tier3_genes.append(gene)

tier3_df = pd.DataFrame({'gene_id': tier3_genes})
tier3_df['tier'] = 3
tier3_df['confidence'] = 'Moderate'
tier3_df['n_methods'] = tier3_df['gene_id'].apply(
    lambda g: sum(1 for gs in gene_sets.values() if g in gs)
)
print(f"  Tier 3 (Moderate): {len(tier3_genes)} genes (2-3 agree)")
print()

# Combine all tiers
tiered_genes = pd.concat([tier1_df, tier2_df, tier3_df], ignore_index=True)
tiered_genes = tiered_genes.sort_values(['tier', 'n_methods'], ascending=[True, False])

print(f"Total consensus genes: {len(tiered_genes)}")
print()

# =============================================================================
# PART 10: Validation Priority List
# =============================================================================

print("PART 10: Creating Validation Priority List")
print("-" * 70)
print()

print("Creating prioritized list for experimental validation...")
print()

# Get Tier 1 genes with direction info
tier1_set = set(tier1_genes)
validation_list = []

for gene in tier1_genes:
    # Get direction from Fisher's result
    gene_data = fisher_result.consensus_genes[
        fisher_result.consensus_genes['gene_id'] == gene
    ]
    
    if len(gene_data) > 0:
        row = gene_data.iloc[0]
        validation_list.append({
            'gene_id': gene,
            'priority': 1,
            'confidence': 'Very High',
            'direction': row.get('direction', 'unknown'),
            'combined_padj': row.get('combined_padj', np.nan),
            'meta_lfc': row.get('meta_lfc', np.nan)
        })

validation_df = pd.DataFrame(validation_list)
validation_df = validation_df.sort_values('combined_padj')

print(f"Validation candidates (Tier 1): {len(validation_df)}")
print()

# Show top 20
print("Top 20 validation candidates:")
print()
print(validation_df.head(20)[['gene_id', 'direction', 
                               'combined_padj']].to_string(index=False))
print()

# =============================================================================
# PART 11: Save Results
# =============================================================================

print("PART 11: Saving Results")
print("-" * 70)
print()

output_dir = Path('output/module9_ensemble/')
output_dir.mkdir(parents=True, exist_ok=True)

# Save comparison
comparison_df.to_csv(output_dir / 'method_comparison.csv', index=False)
print("✓ Saved: method_comparison.csv")

# Save tiered gene list
tiered_genes.to_csv(output_dir / 'tiered_consensus_genes.csv', index=False)
print("✓ Saved: tiered_consensus_genes.csv")

# Save validation list
validation_df.to_csv(output_dir / 'validation_priority_list.csv', index=False)
print("✓ Saved: validation_priority_list.csv")

# Save core consensus
core_df = pd.DataFrame({'gene_id': list(core_consensus)})
core_df.to_csv(output_dir / 'core_consensus_all_methods.csv', index=False)
print("✓ Saved: core_consensus_all_methods.csv")

# Save overlap statistics
overlap_stats = {
    'core_consensus': len(core_consensus),
    'total_unique_genes': len(set.union(*gene_sets.values())),
    'method_counts': {name: len(genes) for name, genes in gene_sets.items()},
    'tier_counts': {
        'tier1': len(tier1_genes),
        'tier2': len(tier2_genes),
        'tier3': len(tier3_genes)
    }
}

import json
with open(output_dir / 'overlap_statistics.json', 'w') as f:
    json.dump(overlap_stats, f, indent=2)
print("✓ Saved: overlap_statistics.json")

print()

# =============================================================================
# PART 12: Recommendations
# =============================================================================

print("PART 12: Recommendations")
print("-" * 70)
print()

# Determine best method
max_genes_method = comparison_df.loc[comparison_df['Genes'].idxmax(), 'Method']
min_genes_method = comparison_df.loc[comparison_df['Genes'].idxmin(), 'Method']

print("Based on your data:")
print()
print("1. FOR MAXIMUM SENSITIVITY (exploratory analysis):")
print(f"   → Use {max_genes_method}")
print(f"   → {comparison_df[comparison_df['Method']==max_genes_method]['Genes'].values[0]} genes")
print()

print("2. FOR HIGH CONFIDENCE (validation experiments):")
print(f"   → Use Tier 1 genes (core consensus)")
print(f"   → {len(tier1_genes)} genes with highest confidence")
print()

print("3. FOR BALANCED APPROACH:")
print("   → Use RRA or Majority Voting")
print(f"   → ~{rra_result.n_consensus_genes} genes")
print()

print("4. FOR PUBLICATION:")
print("   → Report Tier 1 (core consensus) as main results")
print("   → Report Tier 1+2 for comprehensive analysis")
print(f"   → Total: {len(tier1_genes) + len(tier2_genes)} genes")
print()

# =============================================================================
# COMPLETION
# =============================================================================

print("=" * 70)
print("MODULE 9 ENSEMBLE ANALYSIS COMPLETE!")
print("=" * 70)
print()
print("Summary:")
print()
print(f"Analyzed {len(de_results)} DE methods:")
print("  • DESeq2")
print("  • edgeR")
print("  • limma")
print()
print(f"Applied 6 ensemble strategies:")
print("  • Fisher's method: {fisher_result.n_consensus_genes} genes")
print("  • Brown's method: {brown_result.n_consensus_genes} genes")
print("  • RRA: {rra_result.n_consensus_genes} genes")
print("  • Voting (≥2): {voting_majority.n_consensus_genes} genes")
print("  • Voting (≥3): {voting_unanimous.n_consensus_genes} genes")
print("  • Weighted: {weighted_result.n_consensus_genes} genes")
print()
print(f"Created tiered consensus:")
print(f"  • Tier 1 (Very High): {len(tier1_genes)} genes")
print(f"  • Tier 2 (High): {len(tier2_genes)} genes")
print(f"  • Tier 3 (Moderate): {len(tier3_genes)} genes")
print(f"  • Total: {len(tiered_genes)} genes")
print()
print("All results saved to: output/module9_ensemble/")
print()
print("Key files:")
print("  • tiered_consensus_genes.csv - Main gene list")
print("  • validation_priority_list.csv - For qPCR/validation")
print("  • core_consensus_all_methods.csv - Highest confidence")
print("  • method_comparison.csv - Method comparison")
print("  • Individual method directories with detailed results")
print()
print("Next steps:")
print("  1. Use Tier 1 genes for validation experiments")
print("  2. Use Tier 1+2 for pathway enrichment")
print("  3. Use Fisher's results for exploratory analysis")
print("  4. Compare with known biology/literature")
print()
print("=" * 70)

# =============================================================================
# USAGE TIPS
# =============================================================================

print()
print("USAGE TIPS:")
print("-" * 70)
print()

print("Choosing the right method:")
print()
print("  Fisher's/Brown's → Maximum gene discovery")
print("  RRA → Robust ranking, balanced")
print("  Voting (≥3) → Highest confidence, few genes")
print("  Voting (≥2) → Balanced confidence/coverage")
print("  Weighted → Quality-based (requires Module 8)")
print()

print("Parameter tuning:")
print()
print("  # More strict threshold")
print("  result = ensemble_fisher(")
print("      de_results=de_results,")
print("      significance_threshold=0.01  # Instead of 0.05")
print("  )")
print()
print("  # Require higher fold change")
print("  result = ensemble_voting(")
print("      de_results=de_results,")
print("      filters={'padj': 0.05, 'lfc': 1.5}")
print("  )")
print()
print("  # Relax direction consistency")
print("  result = ensemble_fisher(")
print("      de_results=de_results,")
print("      direction_threshold=0.67  # 2/3 majority instead of 1.0")
print("  )")
print()

print("=" * 70)
print("For more information, see: docs/MODULE9_ENSEMBLE_ANALYSIS.md")
print("=" * 70)
