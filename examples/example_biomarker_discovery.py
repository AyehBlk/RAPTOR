#!/usr/bin/env python3
"""
RAPTOR Module 10: Biomarker Discovery — Complete Example Script

Demonstrates ALL features of RAPTOR's biomarker discovery module
using synthetic RNA-seq data.

EXAMPLES COVERED:
    1.  Basic discovery (minimal call, auto panel size)
    2.  Multi-method feature selection with consensus ranking
    3.  Targeted panel size discovery
    4.  Upstream DE result integration (Module 7/9)
    5.  Independent cohort validation
    6.  Component-level usage (fine-grained control)
    7.  Survival biomarker discovery (requires lifelines)
    8.  Dependency status check
    --- Gap-filling examples ---
    9.  Biological annotation (gene info, pathways, literature, PPI, report)
    10. Stability selection
    11. LOOCV with small dataset (< 20 samples)
    12. Save/load BiomarkerResult roundtrip
    13. DataFrame input (not CSV paths)
    14. Feature importance from trained models
    15. Optional methods: Boruta, mRMR, SHAP (if installed)

Run:
    python example_biomarker_discovery.py

Output:
    results/example_*/  — per-example output directories

Author: Ayeh Bolouki
"""

import numpy as np
import pandas as pd
from pathlib import Path
import warnings
warnings.filterwarnings('ignore', category=FutureWarning)


def generate_test_data(output_dir: str = "test_data"):
    """Create synthetic RNA-seq data with known DE genes."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    np.random.seed(42)
    n_genes, n_samples, n_de = 200, 40, 20

    gene_ids = [f"GENE_{i:04d}" for i in range(n_genes)]
    sample_ids = [f"sample_{i:02d}" for i in range(n_samples)]

    counts = np.random.negative_binomial(n=5, p=0.01, size=(n_genes, n_samples))
    for i in range(n_de):
        fold = np.random.uniform(2.0, 5.0)
        counts[i, 20:] = (counts[i, 20:] * fold).astype(int)

    pd.DataFrame(counts, index=gene_ids, columns=sample_ids).to_csv(output_dir / "counts.csv")
    pd.DataFrame({'sample_id': sample_ids, 'condition': ['control']*20 + ['treatment']*20,
                   'batch': [1]*10+[2]*10+[1]*10+[2]*10}).to_csv(output_dir / "metadata.csv", index=False)

    np.random.seed(999)
    val_counts = np.random.negative_binomial(5, 0.01, size=(n_genes, 30))
    val_samples = [f"val_{i:02d}" for i in range(30)]
    for i in range(n_de):
        val_counts[i, 15:] = (val_counts[i, 15:] * 3).astype(int)
    pd.DataFrame(val_counts, index=gene_ids, columns=val_samples).to_csv(output_dir / "validation_counts.csv")
    pd.DataFrame({'sample_id': val_samples, 'condition': ['control']*15+['treatment']*15}).to_csv(
        output_dir / "validation_metadata.csv", index=False)

    np.random.seed(77)
    small_samples = [f"sm_{i}" for i in range(10)]
    small_counts = np.random.negative_binomial(5, 0.01, size=(50, 10))
    for i in range(5):
        small_counts[i, 5:] = (small_counts[i, 5:] * 4).astype(int)
    pd.DataFrame(small_counts, index=[f"G{i}" for i in range(50)], columns=small_samples).to_csv(
        output_dir / "small_counts.csv")
    pd.DataFrame({'sample_id': small_samples, 'condition': ['control']*5+['treatment']*5}).to_csv(
        output_dir / "small_metadata.csv", index=False)

    np.random.seed(123)
    pd.DataFrame({'sample_id': sample_ids,
                   'os_time': np.random.exponential(500, n_samples).astype(int) + 30,
                   'os_event': np.random.binomial(1, 0.6, n_samples)}).to_csv(
        output_dir / "clinical.csv", index=False)

    print(f"Test data generated in {output_dir}/")
    print(f"  counts.csv:              {n_genes} genes x {n_samples} samples")
    print(f"  metadata.csv:            {n_samples} samples, 2 groups")
    print(f"  validation_counts.csv:   {n_genes} genes x 30 samples")
    print(f"  small_counts.csv:        50 genes x 10 samples")
    print(f"  clinical.csv:            {n_samples} patients, survival data")
    print(f"  True DE genes: GENE_0000 .. GENE_0019\n")
    return output_dir


# ── EXAMPLE 1: Basic discovery ──────────────────────────────────────────────

def example_1_basic_discovery(data_dir):
    """Minimal call — auto methods, auto panel size."""
    print("=" * 70)
    print("EXAMPLE 1: Basic Biomarker Discovery")
    print("=" * 70)
    from raptor.biomarker_discovery import discover_biomarkers

    result = discover_biomarkers(
        counts=str(data_dir / "counts.csv"),
        metadata=str(data_dir / "metadata.csv"),
        group_column='condition', annotate=False,
        output_dir='results/example_01_basic', verbose=True)

    best = result.classification_results[result.best_classifier]
    print(f"\nPanel: {result.panel}")
    print(f"Best: {result.best_classifier}, AUC={best.auc:.3f}, F1={best.f1:.3f}\n")
    return result


# ── EXAMPLE 2: Multi-method consensus ───────────────────────────────────────

def example_2_multi_method(data_dir):
    """Multiple methods + consensus ranking."""
    print("=" * 70)
    print("EXAMPLE 2: Multi-Method Feature Selection (target=10)")
    print("=" * 70)
    from raptor.biomarker_discovery import discover_biomarkers

    result = discover_biomarkers(
        counts=str(data_dir / "counts.csv"),
        metadata=str(data_dir / "metadata.csv"),
        group_column='condition', methods=['elastic_net', 'rfe'],
        target_panel_size=10, annotate=False,
        output_dir='results/example_02_multi', verbose=True)

    for m, s in result.selection_results.items():
        print(f"  {m}: {s.n_selected} selected")
    print(f"\nTop 10 consensus:")
    print(result.ranked_genes.head(10)[['consensus_rank','consensus_score','n_methods_selected']])

    de = {f"GENE_{i:04d}" for i in range(20)}
    print(f"\nTrue DE in panel: {len(set(result.panel) & de)}/{result.panel_size}\n")
    return result


# ── EXAMPLE 3: Target panel size ────────────────────────────────────────────

def example_3_target_panel(data_dir):
    """Exact panel size."""
    print("=" * 70)
    print("EXAMPLE 3: Targeted Panel Size (5 genes)")
    print("=" * 70)
    from raptor.biomarker_discovery import discover_biomarkers

    result = discover_biomarkers(
        counts=str(data_dir / "counts.csv"),
        metadata=str(data_dir / "metadata.csv"),
        group_column='condition', methods=['elastic_net'],
        target_panel_size=5, annotate=False,
        output_dir='results/example_03_target', verbose=True)

    best = result.classification_results[result.best_classifier]
    print(f"\nPanel (5): {result.panel}, AUC={best.auc:.3f}\n")
    return result


# ── EXAMPLE 4: Upstream DE integration ──────────────────────────────────────

def example_4_de_integration(data_dir):
    """Use Module 7/9 DE result to restrict candidates."""
    print("=" * 70)
    print("EXAMPLE 4: Upstream DE Integration (Module 7/9)")
    print("=" * 70)
    from raptor.biomarker_discovery import discover_biomarkers

    class MockDEResult:
        significant_genes = [f"GENE_{i:04d}" for i in range(20)]
        n_genes = 20

    result = discover_biomarkers(
        counts=str(data_dir / "counts.csv"),
        metadata=str(data_dir / "metadata.csv"),
        group_column='condition', de_result=MockDEResult(),
        methods=['de_filter', 'elastic_net'], target_panel_size=8,
        annotate=False, output_dir='results/example_04_de', verbose=True)

    de_set = set(MockDEResult.significant_genes)
    print(f"\nPanel: {result.panel}")
    print(f"All from DE list: {all(g in de_set for g in result.panel)}\n")
    return result


# ── EXAMPLE 5: Independent validation ───────────────────────────────────────

def example_5_validation(data_dir, discovery_result):
    """Validate panel on independent cohort."""
    print("=" * 70)
    print("EXAMPLE 5: Independent Cohort Validation")
    print("=" * 70)
    from raptor.biomarker_discovery import validate_biomarkers

    val = validate_biomarkers(
        panel_genes=discovery_result.panel,
        counts=str(data_dir / "validation_counts.csv"),
        metadata=str(data_dir / "validation_metadata.csv"),
        group_column='condition', n_folds=3, verbose=True)

    for name, res in val.items():
        print(f"  {name}: AUC={res.auc:.3f}, F1={res.f1:.3f}")
    best = max(val, key=lambda k: val[k].auc)
    print(f"\nBest validation: {best} (AUC={val[best].auc:.3f})\n")


# ── EXAMPLE 6: Component-level usage ────────────────────────────────────────

def example_6_components(data_dir):
    """Use individual classes for fine-grained control."""
    print("=" * 70)
    print("EXAMPLE 6: Component-Level Usage")
    print("=" * 70)
    from raptor.biomarker_discovery import (
        _prepare_expression_data, FeatureSelector, ClassifierEvaluator, PanelOptimizer)

    counts = pd.read_csv(data_dir / "counts.csv", index_col=0)
    metadata = pd.read_csv(data_dir / "metadata.csv")
    X, y, gene_list = _prepare_expression_data(counts, metadata, 'condition')
    print(f"Prepared: {X.shape[0]} samples, {X.shape[1]} genes")

    selector = FeatureSelector(random_state=42, verbose=False)
    enet = selector.select_lasso(X, y, l1_ratio=0.5, label='elastic_net')
    rfe = selector.select_rfe(X, y, n_features=20)
    print(f"Elastic net: {enet.n_selected}, RFE: {rfe.n_selected}")

    ranking = selector.consensus_ranking(gene_list)
    top_20 = ranking.head(20).index.tolist()

    optimizer = PanelOptimizer(random_state=42, verbose=False)
    panel = optimizer.forward_selection(X, y, ranked_genes=top_20,
                                         min_panel=3, max_panel=15, step=1, n_cv=3)
    print(f"Optimal: {panel.optimal_size} genes, AUC={panel.optimal_auc:.3f}")

    evaluator = ClassifierEvaluator(random_state=42, verbose=False)
    clf = evaluator.evaluate_nested_cv(X[panel.optimal_panel], y, n_outer=3)
    for name, res in clf.items():
        print(f"  {name}: AUC={res.auc:.3f}")
    print()
    return X, y, gene_list


# ── EXAMPLE 7: Survival biomarkers ──────────────────────────────────────────

def example_7_survival(data_dir):
    """Cox regression survival biomarkers (requires lifelines)."""
    print("=" * 70)
    print("EXAMPLE 7: Survival Biomarker Discovery")
    print("=" * 70)
    try:
        from raptor.biomarker_discovery import discover_survival_biomarkers
        result = discover_survival_biomarkers(
            counts=str(data_dir / "counts.csv"),
            clinical=str(data_dir / "clinical.csv"),
            time_column='os_time', event_column='os_event',
            fdr_threshold=0.5, output_dir='results/example_07_survival', verbose=True)
        print(f"\nPrognostic genes: {len(result.significant_genes)}, C-index: {result.c_index:.3f}")
    except ImportError:
        print("  SKIPPED — pip install lifelines")
    except Exception as e:
        print(f"  Note: {e} (expected with synthetic data)")
    print()


# ── EXAMPLE 8: Dependency check ─────────────────────────────────────────────

def example_8_check_deps():
    """Check optional dependency availability."""
    print("=" * 70)
    print("EXAMPLE 8: Dependency Status")
    print("=" * 70)
    from raptor.biomarker_discovery import get_dependencies_status
    for pkg, ok in get_dependencies_status().items():
        print(f"  {'✓' if ok else '✗'} {pkg}")
    print()


# ── EXAMPLE 9: Biological annotation ────────────────────────────────────────

def example_9_annotation(data_dir):
    """Full annotation pipeline: gene info, pathways, literature, PPI, report."""
    print("=" * 70)
    print("EXAMPLE 9: Biological Annotation & Report Generation")
    print("=" * 70)
    from raptor.biomarker_discovery import (
        BiologicalAnnotator, AnnotationResult, ClassificationResult, PanelOptimizationResult)

    output_dir = Path('results/example_09_annotation')
    output_dir.mkdir(parents=True, exist_ok=True)
    panel_genes = ['TP53', 'BRCA1', 'EGFR', 'KRAS', 'MYC']
    annotator = BiologicalAnnotator(species='human', verbose=True)

    # 9A: Gene annotation
    print("\n--- 9A: Gene Annotation (MyGene.info) ---")
    gene_info = pd.DataFrame()
    try:
        gene_info = annotator.annotate_genes(panel_genes)
        if not gene_info.empty:
            for gid, row in gene_info.iterrows():
                print(f"  {gid}: {row.get('name', 'N/A')[:60]}")
    except Exception as e:
        print(f"  Failed: {e}")

    # 9B: Pathway enrichment
    print("\n--- 9B: Pathway Enrichment ---")
    enrichment = pd.DataFrame()
    try:
        enrichment = annotator.pathway_enrichment(panel_genes)
        if not enrichment.empty:
            print(f"  {len(enrichment)} pathways found")
            print(enrichment.head(5)[['pathway', 'p_value']].to_string())
    except Exception as e:
        print(f"  Failed: {e}")

    # 9C: Literature search
    print("\n--- 9C: Literature Mining (Europe PMC) ---")
    lit = pd.DataFrame()
    try:
        lit = annotator.literature_search(panel_genes, disease_term='cancer')
        if not lit.empty:
            for _, row in lit.iterrows():
                print(f"  {row['gene_id']}: {row['n_publications']} publications")
    except Exception as e:
        print(f"  Failed: {e}")

    # 9D: PPI network
    print("\n--- 9D: PPI Network (STRING) ---")
    ppi = None
    try:
        ppi = annotator.ppi_network(panel_genes)
        if ppi:
            print(f"  Nodes: {ppi['n_nodes']}, Edges: {ppi['n_edges']}")
            if ppi.get('network_url'):
                print(f"  URL: {ppi['network_url']}")
    except Exception as e:
        print(f"  Failed: {e}")

    # 9E: Full pipeline
    print("\n--- 9E: Full annotate_panel() ---")
    ann_result = None
    try:
        ann_result = annotator.annotate_panel(panel_genes, disease_term='cancer')
        print(f"  Annotated: {ann_result.n_annotated}, Pathways: {ann_result.n_enriched_pathways}")
        ann_result.save(output_dir)
    except Exception as e:
        print(f"  Failed: {e}")
        ann_result = AnnotationResult(gene_annotations=gene_info,
                                       pathway_enrichment=enrichment,
                                       literature_hits=lit, ppi_network=ppi)

    # 9F: Report generation
    print("\n--- 9F: Markdown Report ---")
    try:
        mock_clf = {'random_forest': ClassificationResult(
            model_name='random_forest', auc=0.95, f1=0.92, sensitivity=0.93, specificity=0.90)}
        mock_opt = PanelOptimizationResult(
            optimal_panel=panel_genes, optimal_size=5, optimal_auc=0.95,
            panel_curve=pd.DataFrame({'panel_size': [3,4,5], 'auc_mean': [0.88,0.92,0.95]}))

        report = annotator.generate_report(
            panel_genes=panel_genes, annotation_result=ann_result,
            classification_results=mock_clf, panel_optimization=mock_opt,
            output_path=output_dir / "biomarker_report.md")
        print(f"  Report: {len(report)} chars -> {output_dir}/biomarker_report.md")
        for line in report.split('\n')[:8]:
            print(f"    {line}")
    except Exception as e:
        print(f"  Failed: {e}")
    print()


# ── EXAMPLE 10: Stability selection ─────────────────────────────────────────

def example_10_stability(data_dir):
    """Bootstrap stability selection — which genes are consistently chosen?"""
    print("=" * 70)
    print("EXAMPLE 10: Stability Selection")
    print("=" * 70)
    from raptor.biomarker_discovery import _prepare_expression_data, PanelOptimizer

    counts = pd.read_csv(data_dir / "counts.csv", index_col=0)
    metadata = pd.read_csv(data_dir / "metadata.csv")
    X, y, _ = _prepare_expression_data(counts, metadata, 'condition')

    optimizer = PanelOptimizer(random_state=42, verbose=True)
    stab = optimizer.stability_selection(X, y, n_bootstrap=20, threshold=0.3)

    n_stable = stab['is_stable'].sum()
    print(f"\nStable genes (freq >= 0.3): {n_stable}")
    print(f"\nTop 10:")
    for gid, row in stab.nlargest(10, 'selection_frequency').iterrows():
        de = "  (DE)" if int(gid.split('_')[1]) < 20 else ""
        print(f"  {gid}: freq={row['selection_frequency']:.2f}{de}")
    print()


# ── EXAMPLE 11: LOOCV small dataset ─────────────────────────────────────────

def example_11_loocv_small(data_dir):
    """Leave-one-out CV for small datasets (< 20 samples)."""
    print("=" * 70)
    print("EXAMPLE 11: LOOCV with Small Dataset (10 samples)")
    print("=" * 70)
    from raptor.biomarker_discovery import (
        _prepare_expression_data, FeatureSelector, ClassifierEvaluator, discover_biomarkers)

    counts = pd.read_csv(data_dir / "small_counts.csv", index_col=0)
    metadata = pd.read_csv(data_dir / "small_metadata.csv")
    X, y, _ = _prepare_expression_data(counts, metadata, 'condition')
    print(f"Small dataset: {X.shape[0]} samples, {X.shape[1]} genes")

    # Explicit LOOCV
    selector = FeatureSelector(random_state=42, verbose=False)
    enet = selector.select_lasso(X, y, l1_ratio=0.5, label='elastic_net')
    top5 = enet.gene_scores.nlargest(5, 'score').index.tolist()

    evaluator = ClassifierEvaluator(random_state=42, verbose=False)
    results = evaluator.evaluate_loocv(X[top5], y)
    print(f"\nExplicit LOOCV:")
    for name, res in results.items():
        print(f"  {name}: AUC={res.auc:.3f}, F1={res.f1:.3f}")

    # Auto-LOOCV via discover_biomarkers
    result = discover_biomarkers(
        counts=str(data_dir / "small_counts.csv"),
        metadata=str(data_dir / "small_metadata.csv"),
        group_column='condition', methods=['elastic_net'], target_panel_size=3,
        annotate=False, output_dir='results/example_11_loocv', verbose=False)
    best = result.classification_results[result.best_classifier]
    print(f"\nAuto-adapt: panel={result.panel}, AUC={best.auc:.3f}\n")


# ── EXAMPLE 12: Save/Load roundtrip ─────────────────────────────────────────

def example_12_save_load(data_dir):
    """Pickle persistence and reload."""
    print("=" * 70)
    print("EXAMPLE 12: Save/Load BiomarkerResult Roundtrip")
    print("=" * 70)
    from raptor.biomarker_discovery import discover_biomarkers, BiomarkerResult

    out = Path('results/example_12_saveload')
    original = discover_biomarkers(
        counts=str(data_dir / "counts.csv"),
        metadata=str(data_dir / "metadata.csv"),
        group_column='condition', methods=['elastic_net'], target_panel_size=5,
        annotate=False, output_dir=str(out), verbose=False)

    best_auc = original.classification_results[original.best_classifier].auc
    print(f"Original: panel={original.panel}, AUC={best_auc:.3f}")

    print(f"\nSaved files:")
    for f in sorted(out.iterdir()):
        print(f"  {f.name:40s} ({f.stat().st_size:,} bytes)")

    loaded = BiomarkerResult.load(out)
    assert loaded.panel == original.panel
    assert loaded.best_classifier == original.best_classifier
    print(f"\n✓ Roundtrip verified")

    print(f"\nLoaded summary:")
    print(loaded.summary())
    print()


# ── EXAMPLE 13: DataFrame input ─────────────────────────────────────────────

def example_13_dataframe_input(data_dir):
    """Pass DataFrames instead of file paths."""
    print("=" * 70)
    print("EXAMPLE 13: DataFrame Input (In-Memory)")
    print("=" * 70)
    from raptor.biomarker_discovery import discover_biomarkers

    counts_df = pd.read_csv(data_dir / "counts.csv", index_col=0)
    metadata_df = pd.read_csv(data_dir / "metadata.csv")
    print(f"In-memory: {counts_df.shape[0]} genes x {counts_df.shape[1]} samples")

    result = discover_biomarkers(
        counts=counts_df, metadata=metadata_df,  # DataFrames, not paths
        group_column='condition', methods=['elastic_net'], target_panel_size=5,
        annotate=False, output_dir='results/example_13_dataframe', verbose=False)

    best = result.classification_results[result.best_classifier]
    print(f"Panel: {result.panel}, AUC={best.auc:.3f}\n")


# ── EXAMPLE 14: Feature importance ──────────────────────────────────────────

def example_14_feature_importance(data_dir):
    """Access feature importance & trained models."""
    print("=" * 70)
    print("EXAMPLE 14: Feature Importance from Trained Models")
    print("=" * 70)
    from raptor.biomarker_discovery import _prepare_expression_data, ClassifierEvaluator

    counts = pd.read_csv(data_dir / "counts.csv", index_col=0)
    metadata = pd.read_csv(data_dir / "metadata.csv")
    X, y, gene_list = _prepare_expression_data(counts, metadata, 'condition')

    panel = [g for g in gene_list if int(g.split('_')[1]) < 20][:10]
    evaluator = ClassifierEvaluator(random_state=42, verbose=False)
    results = evaluator.evaluate_nested_cv(X[panel], y, n_outer=3,
                                            classifiers=['random_forest', 'logistic_regression'])

    # Feature importance (Random Forest)
    rf = results['random_forest']
    print("Random Forest feature importance:")
    if rf.feature_importance is not None:
        imp = rf.feature_importance.sort_values('importance', ascending=False)
        for gid, row in imp.iterrows():
            bar = "█" * int(row['importance'] * 50)
            print(f"  {gid}: {row['importance']:.4f} {bar}")

    # Trained model access
    print(f"\nTrained model: {type(rf.trained_model).__name__}, "
          f"n_estimators={rf.trained_model.n_estimators}")

    # Per-fold metrics
    print(f"\nPer-fold metrics:")
    for i, fold in enumerate(rf.metrics_per_fold):
        print(f"  Fold {i+1}: AUC={fold.get('auc', 0):.3f}, F1={fold.get('f1', 0):.3f}")
    print()


# ── EXAMPLE 15: Optional methods ────────────────────────────────────────────

def example_15_optional_methods(data_dir):
    """Boruta, mRMR, SHAP — run if installed, skip if not."""
    print("=" * 70)
    print("EXAMPLE 15: Optional Feature Selection Methods")
    print("=" * 70)
    from raptor.biomarker_discovery import (
        _prepare_expression_data, FeatureSelector,
        _BORUTA_AVAILABLE, _MRMR_AVAILABLE, _SHAP_AVAILABLE)

    counts = pd.read_csv(data_dir / "counts.csv", index_col=0)
    metadata = pd.read_csv(data_dir / "metadata.csv")
    X, y, gene_list = _prepare_expression_data(counts, metadata, 'condition')
    selector = FeatureSelector(random_state=42, verbose=True)

    print(f"\n--- Boruta (available: {_BORUTA_AVAILABLE}) ---")
    if _BORUTA_AVAILABLE:
        r = selector.select_boruta(X, y, max_iter=30)
        print(f"  Selected: {r.n_selected}, top: {r.selected_genes[:5]}")
    else:
        print("  SKIPPED — pip install boruta")

    print(f"\n--- mRMR (available: {_MRMR_AVAILABLE}) ---")
    if _MRMR_AVAILABLE:
        r = selector.select_mrmr(X, y, n_features=15)
        print(f"  Selected: {r.n_selected}, top: {r.selected_genes[:5]}")
    else:
        print("  SKIPPED — pip install mrmr-selection")

    print(f"\n--- SHAP (available: {_SHAP_AVAILABLE}) ---")
    if _SHAP_AVAILABLE:
        r = selector.select_shap(X, y, n_features=15)
        print(f"  Selected: {r.n_selected}, top: {r.selected_genes[:5]}")
    else:
        print("  SKIPPED — pip install shap")

    print(f"\nMethods run: {len(selector.results)}")
    for m, r in selector.results.items():
        print(f"  {m}: {r.n_selected} selected")

    if len(selector.results) > 1:
        ranking = selector.consensus_ranking(gene_list)
        print(f"\nConsensus top 5:")
        for gid, row in ranking.head(5).iterrows():
            print(f"  {gid}: rank={row['consensus_rank']}, methods={row['n_methods_selected']}")
    print()


# ── MAIN ────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    print("\n🦖 RAPTOR Module 10: Biomarker Discovery — Complete Examples")
    print("=" * 70 + "\n")

    data_dir = generate_test_data("test_data")
    example_8_check_deps()

    r1 = example_1_basic_discovery(data_dir)
    r2 = example_2_multi_method(data_dir)
    r3 = example_3_target_panel(data_dir)
    r4 = example_4_de_integration(data_dir)
    example_5_validation(data_dir, r2)
    X, y, gl = example_6_components(data_dir)
    example_7_survival(data_dir)
    example_9_annotation(data_dir)
    example_10_stability(data_dir)
    example_11_loocv_small(data_dir)
    example_12_save_load(data_dir)
    example_13_dataframe_input(data_dir)
    example_14_feature_importance(data_dir)
    example_15_optional_methods(data_dir)

    print("\n" + "=" * 70)
    print("ALL 15 EXAMPLES COMPLETED!")
    print("=" * 70)
    print("\nOutput directories:")
    for d in sorted(Path("results").iterdir()):
        if d.is_dir() and d.name.startswith("example_"):
            print(f"  {d.name}/  ({len(list(d.rglob('*')))} files)")
    print("\nFeatures covered:")
    for f in [
        "10A: Elastic Net, LASSO, RFE, DE filter, Boruta*, mRMR*, SHAP*",
        "10A: Consensus ranking across methods",
        "10B: Nested CV (LogReg, RF, SVM, XGBoost*)",
        "10B: LOOCV for small datasets (explicit + auto)",
        "10B: Feature importance & trained model access",
        "10B: Per-fold metrics access",
        "10C: Forward panel selection, auto & target size",
        "10C: Stability selection (bootstrap)",
        "10D: Cox univariate screen & CoxNet panel*",
        "Annotation: Gene info (MyGene), pathways, literature, PPI",
        "Annotation: Full annotate_panel() pipeline",
        "Annotation: Markdown report generation",
        "I/O: CSV path input & DataFrame input",
        "I/O: Save/load BiomarkerResult roundtrip",
        "Workflow: Discovery -> independent validation",
        "Workflow: Upstream DE integration (M7/M9)",
        "* = requires optional dependency",
    ]:
        print(f"  ✓ {f}")
    print()
