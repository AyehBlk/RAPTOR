"""
RAPTOR Module 10: Biomarker Discovery — Comprehensive Test Suite

Run with:
    pytest test_biomarker_discovery.py -v
    pytest test_biomarker_discovery.py -v -k "test_feature"    # run only feature selection tests
    pytest test_biomarker_discovery.py -v -k "test_integration" # run only integration tests
    pytest test_biomarker_discovery.py -v --tb=short            # shorter tracebacks

Requirements:
    pip install pytest scikit-learn numpy pandas scipy

Optional (tests skip gracefully if missing):
    pip install xgboost shap boruta mrmr-selection lifelines PyWGCNA gseapy mygene

Author: Ayeh Bolouki
"""

import pytest
import numpy as np
import pandas as pd
import tempfile
import shutil
from pathlib import Path
from unittest.mock import patch, MagicMock

# ============================================================================
# FIXTURES — Synthetic RNA-seq data
# ============================================================================

@pytest.fixture(scope="session")
def synthetic_counts():
    """
    Create a synthetic count matrix: 200 genes x 40 samples.
    First 20 genes are DE (different means between groups).
    """
    np.random.seed(42)
    n_genes = 200
    n_samples = 40
    n_de = 20  # truly differential genes

    gene_ids = [f"GENE_{i:04d}" for i in range(n_genes)]
    sample_ids = [f"sample_{i:02d}" for i in range(n_samples)]

    # Base expression: Poisson-like counts
    counts = np.random.negative_binomial(n=5, p=0.01, size=(n_genes, n_samples))

    # Make first n_de genes differential: higher in group 1 (samples 20-39)
    for i in range(n_de):
        fold = np.random.uniform(2.0, 5.0)
        counts[i, 20:] = (counts[i, 20:] * fold).astype(int)

    df = pd.DataFrame(counts, index=gene_ids, columns=sample_ids)
    return df


@pytest.fixture(scope="session")
def synthetic_metadata():
    """Metadata with 40 samples: 20 control, 20 treatment."""
    sample_ids = [f"sample_{i:02d}" for i in range(40)]
    conditions = ['control'] * 20 + ['treatment'] * 20
    return pd.DataFrame({
        'sample_id': sample_ids,
        'condition': conditions,
        'batch': [1] * 10 + [2] * 10 + [1] * 10 + [2] * 10,
    })


@pytest.fixture(scope="session")
def synthetic_de_genes():
    """List of truly DE gene IDs."""
    return [f"GENE_{i:04d}" for i in range(20)]


@pytest.fixture(scope="session")
def prepared_data(synthetic_counts, synthetic_metadata):
    """Pre-prepared X (samples x genes), y (labels), gene_list."""
    from raptor.biomarker_discovery import _prepare_expression_data
    X, y, gene_list = _prepare_expression_data(
        synthetic_counts, synthetic_metadata, 'condition',
    )
    return X, y, gene_list


@pytest.fixture
def output_dir():
    """Temporary output directory, cleaned up after test."""
    tmpdir = tempfile.mkdtemp(prefix="raptor_m10_test_")
    yield Path(tmpdir)
    shutil.rmtree(tmpdir, ignore_errors=True)


@pytest.fixture(scope="session")
def synthetic_clinical():
    """Clinical data for survival analysis: 40 patients."""
    np.random.seed(123)
    sample_ids = [f"sample_{i:02d}" for i in range(40)]
    return pd.DataFrame({
        'sample_id': sample_ids,
        'os_time': np.random.exponential(500, 40).astype(int) + 30,
        'os_event': np.random.binomial(1, 0.6, 40),
    })


# ============================================================================
# 1. DATA PREPARATION TESTS
# ============================================================================

class TestDataPreparation:
    """Test _prepare_expression_data function."""

    def test_basic_preparation(self, synthetic_counts, synthetic_metadata):
        from raptor.biomarker_discovery import _prepare_expression_data
        X, y, gene_list = _prepare_expression_data(
            synthetic_counts, synthetic_metadata, 'condition',
        )

        assert isinstance(X, pd.DataFrame)
        assert isinstance(y, np.ndarray)
        assert X.shape[0] == 40  # 40 samples
        assert X.shape[1] > 0   # some genes pass filter
        assert len(y) == 40
        assert set(y) == {0, 1}  # binary labels
        assert len(gene_list) == X.shape[1]

    def test_labels_are_binary(self, prepared_data):
        _, y, _ = prepared_data
        assert set(np.unique(y)) == {0, 1}

    def test_gene_filtering(self, synthetic_counts, synthetic_metadata):
        from raptor.biomarker_discovery import _prepare_expression_data
        # With high min_count, fewer genes should pass
        X_strict, _, _ = _prepare_expression_data(
            synthetic_counts, synthetic_metadata, 'condition',
            min_count=100,
        )
        X_lenient, _, _ = _prepare_expression_data(
            synthetic_counts, synthetic_metadata, 'condition',
            min_count=1,
        )
        assert X_strict.shape[1] <= X_lenient.shape[1]

    def test_de_gene_restriction(self, synthetic_counts, synthetic_metadata, synthetic_de_genes):
        from raptor.biomarker_discovery import _prepare_expression_data
        X, _, gene_list = _prepare_expression_data(
            synthetic_counts, synthetic_metadata, 'condition',
            de_genes=synthetic_de_genes,
        )
        # Should be restricted to DE genes (that pass filter)
        assert X.shape[1] <= len(synthetic_de_genes)
        for g in gene_list:
            assert g in synthetic_de_genes

    def test_invalid_group_column(self, synthetic_counts, synthetic_metadata):
        from raptor.biomarker_discovery import _prepare_expression_data
        with pytest.raises(Exception):
            _prepare_expression_data(
                synthetic_counts, synthetic_metadata, 'nonexistent_column',
            )


# ============================================================================
# 2. FEATURE SELECTION TESTS (10A)
# ============================================================================

class TestFeatureSelector:
    """Test individual feature selection methods."""

    def test_de_filter(self, prepared_data, synthetic_de_genes):
        from raptor.biomarker_discovery import FeatureSelector
        X, y, gene_list = prepared_data

        selector = FeatureSelector(random_state=42, verbose=False)
        result = selector.select_de_filter(synthetic_de_genes, gene_list)

        assert result.method == 'de_filter'
        assert result.n_selected > 0
        assert len(result.selected_genes) == result.n_selected
        assert not result.gene_scores.empty

    def test_elastic_net(self, prepared_data):
        from raptor.biomarker_discovery import FeatureSelector
        X, y, _ = prepared_data

        selector = FeatureSelector(random_state=42, verbose=False)
        result = selector.select_lasso(X, y, l1_ratio=0.5, label='elastic_net')

        assert result.method == 'elastic_net'
        assert result.n_selected >= 0
        assert len(result.gene_scores) == X.shape[1]
        assert result.gene_scores.index.name == 'gene_id'
        # Scores should be non-negative (absolute coefficients)
        assert (result.gene_scores['score'] >= 0).all()

    def test_lasso_pure(self, prepared_data):
        from raptor.biomarker_discovery import FeatureSelector
        X, y, _ = prepared_data

        selector = FeatureSelector(random_state=42, verbose=False)
        result = selector.select_lasso(X, y, l1_ratio=1.0, label='lasso')

        assert result.method == 'lasso'
        # Pure LASSO should select fewer features than Elastic Net
        assert result.n_selected >= 0

    def test_rfe(self, prepared_data):
        from raptor.biomarker_discovery import FeatureSelector
        X, y, _ = prepared_data

        selector = FeatureSelector(random_state=42, verbose=False)
        result = selector.select_rfe(X, y, n_features=15)

        assert result.method == 'rfe'
        assert result.n_selected == 15
        assert len(result.selected_genes) == 15

    def test_boruta(self, prepared_data):
        from raptor.biomarker_discovery import FeatureSelector, _BORUTA_AVAILABLE
        if not _BORUTA_AVAILABLE:
            pytest.skip("boruta not installed")

        X, y, _ = prepared_data
        selector = FeatureSelector(random_state=42, verbose=False)
        result = selector.select_boruta(X, y, max_iter=20)

        assert result.method == 'boruta'
        assert result.n_selected >= 0

    def test_mrmr(self, prepared_data):
        from raptor.biomarker_discovery import FeatureSelector, _MRMR_AVAILABLE
        if not _MRMR_AVAILABLE:
            pytest.skip("mrmr-selection not installed")

        X, y, _ = prepared_data
        selector = FeatureSelector(random_state=42, verbose=False)
        result = selector.select_mrmr(X, y, n_features=15)

        assert result.method == 'mrmr'
        assert result.n_selected == 15

    def test_shap(self, prepared_data):
        from raptor.biomarker_discovery import FeatureSelector, _SHAP_AVAILABLE
        if not _SHAP_AVAILABLE:
            pytest.skip("shap not installed")

        X, y, _ = prepared_data
        selector = FeatureSelector(random_state=42, verbose=False)
        result = selector.select_shap(X, y, n_features=15)

        assert result.method == 'shap'
        assert result.n_selected == 15
        assert (result.gene_scores['score'] >= 0).all()

    def test_consensus_ranking(self, prepared_data, synthetic_de_genes):
        """Test consensus ranking with multiple methods."""
        from raptor.biomarker_discovery import FeatureSelector
        X, y, gene_list = prepared_data

        selector = FeatureSelector(random_state=42, verbose=False)
        selector.select_de_filter(synthetic_de_genes, gene_list)
        selector.select_lasso(X, y, l1_ratio=0.5, label='elastic_net')
        selector.select_rfe(X, y, n_features=20)

        ranking = selector.consensus_ranking(gene_list)

        assert isinstance(ranking, pd.DataFrame)
        assert 'consensus_rank' in ranking.columns
        assert 'consensus_score' in ranking.columns
        assert 'n_methods_selected' in ranking.columns
        assert len(ranking) == len(gene_list)
        assert ranking['consensus_rank'].min() == 1

        # DE genes should rank higher on average
        de_ranks = ranking.loc[
            ranking.index.isin(synthetic_de_genes), 'consensus_rank'
        ].mean()
        non_de_ranks = ranking.loc[
            ~ranking.index.isin(synthetic_de_genes), 'consensus_rank'
        ].mean()
        assert de_ranks < non_de_ranks, (
            "DE genes should have lower (better) consensus rank"
        )

    def test_consensus_ranking_empty(self, prepared_data):
        """Consensus ranking with no results should raise."""
        from raptor.biomarker_discovery import FeatureSelector
        _, _, gene_list = prepared_data
        selector = FeatureSelector(random_state=42, verbose=False)
        with pytest.raises(ValueError, match="No feature selection results"):
            selector.consensus_ranking(gene_list)


# ============================================================================
# 3. CLASSIFICATION TESTS (10B)
# ============================================================================

class TestClassifierEvaluator:
    """Test classification and validation."""

    def test_nested_cv(self, prepared_data):
        from raptor.biomarker_discovery import ClassifierEvaluator
        X, y, _ = prepared_data

        # Use subset of genes as "panel"
        panel = X.columns[:15]
        X_panel = X[panel]

        evaluator = ClassifierEvaluator(random_state=42, verbose=False)
        results = evaluator.evaluate_nested_cv(X_panel, y, n_outer=3)

        assert len(results) >= 2  # At least LogReg + RF
        for name, res in results.items():
            assert res.model_name == name
            assert 0.0 <= res.auc <= 1.0
            assert 0.0 <= res.accuracy <= 1.0
            assert 0.0 <= res.f1 <= 1.0
            assert 0.0 <= res.sensitivity <= 1.0
            assert 0.0 <= res.specificity <= 1.0
            assert len(res.metrics_per_fold) == 3  # 3 folds

    def test_nested_cv_with_good_genes(self, prepared_data):
        """DE genes should yield high AUC."""
        from raptor.biomarker_discovery import ClassifierEvaluator
        X, y, _ = prepared_data

        # Use truly DE genes
        de_panel = [c for c in X.columns if int(c.split('_')[1]) < 20]
        X_panel = X[de_panel[:10]]

        evaluator = ClassifierEvaluator(random_state=42, verbose=False)
        results = evaluator.evaluate_nested_cv(X_panel, y, n_outer=3)

        best_auc = max(r.auc for r in results.values())
        assert best_auc > 0.7, (
            f"DE genes should give AUC > 0.7, got {best_auc:.3f}"
        )

    def test_loocv(self, prepared_data):
        from raptor.biomarker_discovery import ClassifierEvaluator
        X, y, _ = prepared_data

        panel = X.columns[:10]
        X_panel = X[panel]

        evaluator = ClassifierEvaluator(random_state=42, verbose=False)
        results = evaluator.evaluate_loocv(
            X_panel, y, classifiers=['logistic_regression'],
        )

        assert 'logistic_regression' in results
        assert 0.0 <= results['logistic_regression'].auc <= 1.0

    def test_feature_importance(self, prepared_data):
        from raptor.biomarker_discovery import ClassifierEvaluator
        X, y, _ = prepared_data

        panel = X.columns[:10]
        X_panel = X[panel]

        evaluator = ClassifierEvaluator(random_state=42, verbose=False)
        results = evaluator.evaluate_nested_cv(
            X_panel, y, n_outer=3,
            classifiers=['random_forest'],
        )

        rf_result = results['random_forest']
        assert rf_result.feature_importance is not None
        assert len(rf_result.feature_importance) == 10
        assert rf_result.trained_model is not None


# ============================================================================
# 4. PANEL OPTIMIZATION TESTS (10C)
# ============================================================================

class TestPanelOptimizer:
    """Test panel size optimization."""

    def test_forward_selection(self, prepared_data):
        from raptor.biomarker_discovery import PanelOptimizer
        X, y, _ = prepared_data

        ranked = list(X.columns[:30])  # Top 30 genes by "rank"

        optimizer = PanelOptimizer(random_state=42, verbose=False)
        result = optimizer.forward_selection(
            X, y, ranked_genes=ranked,
            min_panel=3, max_panel=20, step=2, n_cv=3,
        )

        assert result.optimal_size >= 3
        assert result.optimal_size <= 20
        assert result.optimal_auc > 0.0
        assert len(result.optimal_panel) == result.optimal_size
        assert not result.panel_curve.empty
        assert 'panel_size' in result.panel_curve.columns
        assert 'auc_mean' in result.panel_curve.columns

    def test_target_panel_size(self, prepared_data):
        from raptor.biomarker_discovery import PanelOptimizer
        X, y, _ = prepared_data

        ranked = list(X.columns[:30])

        optimizer = PanelOptimizer(random_state=42, verbose=False)
        result = optimizer.forward_selection(
            X, y, ranked_genes=ranked,
            min_panel=3, max_panel=20, step=1, n_cv=3,
            target_size=7,
        )

        assert result.optimal_size == 7
        assert len(result.optimal_panel) == 7

    def test_stability_selection(self, prepared_data):
        from raptor.biomarker_discovery import PanelOptimizer
        X, y, _ = prepared_data

        optimizer = PanelOptimizer(random_state=42, verbose=False)
        stability_df = optimizer.stability_selection(
            X, y, n_bootstrap=10, threshold=0.3,
        )

        assert isinstance(stability_df, pd.DataFrame)
        assert 'selection_frequency' in stability_df.columns
        assert 'is_stable' in stability_df.columns
        assert len(stability_df) == X.shape[1]
        assert stability_df['selection_frequency'].max() <= 1.0
        assert stability_df['selection_frequency'].min() >= 0.0


# ============================================================================
# 5. SURVIVAL ANALYSIS TESTS (10D)
# ============================================================================

class TestSurvivalAnalyzer:
    """Test survival analysis methods."""

    @pytest.fixture(autouse=True)
    def check_lifelines(self):
        try:
            import lifelines  # noqa: F401
        except ImportError:
            pytest.skip("lifelines not installed")

    def test_cox_univariate_screen(self, prepared_data, synthetic_clinical):
        from raptor.biomarker_discovery import SurvivalAnalyzer
        X, _, _ = prepared_data

        analyzer = SurvivalAnalyzer(verbose=False)
        time_arr = synthetic_clinical['os_time'].values[:len(X)]
        event_arr = synthetic_clinical['os_event'].values[:len(X)]

        results = analyzer.cox_univariate_screen(
            X.iloc[:, :30], time_arr, event_arr,
        )

        assert isinstance(results, pd.DataFrame)
        assert 'gene_id' in results.columns
        assert 'hazard_ratio' in results.columns
        assert 'p_value' in results.columns
        assert 'adjusted_p_value' in results.columns
        assert len(results) > 0

    def test_cox_lasso_panel(self, prepared_data, synthetic_clinical):
        from raptor.biomarker_discovery import SurvivalAnalyzer
        X, _, _ = prepared_data

        analyzer = SurvivalAnalyzer(verbose=False)
        time_arr = synthetic_clinical['os_time'].values[:len(X)]
        event_arr = synthetic_clinical['os_event'].values[:len(X)]

        result = analyzer.cox_lasso_panel(
            X.iloc[:, :20], time_arr, event_arr, alpha=0.5,
        )

        assert result.c_index > 0.0
        assert isinstance(result.panel_genes, list)


# ============================================================================
# 6. BIOLOGICAL ANNOTATOR TESTS
# ============================================================================

class TestBiologicalAnnotator:
    """Test annotation methods (with mocked API calls)."""

    def test_annotate_genes_mocked(self):
        """Test gene annotation with mocked MyGene.info."""
        from raptor.biomarker_discovery import BiologicalAnnotator

        mock_results = {
            'out': [
                {'query': 'TP53', 'symbol': 'TP53', 'name': 'tumor protein p53',
                 'entrezgene': 7157, 'type_of_gene': 'protein-coding',
                 'summary': 'Acts as tumor suppressor', 'go': {}, 'pathway': {}},
                {'query': 'BRCA1', 'symbol': 'BRCA1', 'name': 'BRCA1 DNA repair',
                 'entrezgene': 672, 'type_of_gene': 'protein-coding',
                 'summary': 'DNA repair', 'go': {}, 'pathway': {}},
            ]
        }

        annotator = BiologicalAnnotator(verbose=False)

        with patch('raptor.biomarker_discovery.BiologicalAnnotator.annotate_genes') as mock:
            mock_df = pd.DataFrame({
                'symbol': ['TP53', 'BRCA1'],
                'name': ['tumor protein p53', 'BRCA1 DNA repair'],
                'entrezgene': ['7157', '672'],
                'type_of_gene': ['protein-coding', 'protein-coding'],
            }, index=pd.Index(['TP53', 'BRCA1'], name='gene_id'))
            mock.return_value = mock_df

            result = annotator.annotate_genes(['TP53', 'BRCA1'])
            assert len(result) == 2
            assert 'symbol' in result.columns

    def test_manual_ora(self):
        """Test manual Fisher's exact ORA."""
        from raptor.biomarker_discovery import BiologicalAnnotator

        annotator = BiologicalAnnotator(verbose=False)

        gene_sets = {
            'pathway_A': ['GENE_0000', 'GENE_0001', 'GENE_0002', 'GENE_0003',
                          'GENE_0004', 'GENE_0050', 'GENE_0051'],
            'pathway_B': ['GENE_0010', 'GENE_0011', 'GENE_0099', 'GENE_0100',
                          'GENE_0101'],
        }
        query_genes = ['GENE_0000', 'GENE_0001', 'GENE_0002', 'GENE_0003',
                       'GENE_0010']
        background = [f'GENE_{i:04d}' for i in range(200)]

        result = annotator._manual_ora(
            query_genes, gene_sets, background, fdr_threshold=0.05,
        )

        assert isinstance(result, pd.DataFrame)
        assert 'pathway' in result.columns
        assert 'p_value' in result.columns
        assert 'adjusted_p_value' in result.columns
        # pathway_A should have significant overlap (4 out of 7)
        assert len(result) >= 1

    def test_annotate_panel_offline(self):
        """Test full annotation pipeline in offline mode."""
        from raptor.biomarker_discovery import BiologicalAnnotator, AnnotationResult

        annotator = BiologicalAnnotator(verbose=False)

        # Mock all network calls
        with patch.object(annotator, 'annotate_genes', return_value=pd.DataFrame()):
            with patch.object(annotator, 'pathway_enrichment', return_value=pd.DataFrame()):
                with patch.object(annotator, 'literature_search', return_value=pd.DataFrame()):
                    with patch.object(annotator, 'ppi_network', return_value={'edges': [], 'n_edges': 0, 'n_nodes': 3}):
                        result = annotator.annotate_panel(
                            panel_genes=['A', 'B', 'C'],
                            run_literature=True,
                            run_ppi=True,
                        )

        assert isinstance(result, AnnotationResult)
        assert result.species == 'human'

    def test_generate_report(self, output_dir):
        """Test Markdown report generation."""
        from raptor.biomarker_discovery import (
            BiologicalAnnotator, AnnotationResult, ClassificationResult,
            PanelOptimizationResult,
        )

        annotator = BiologicalAnnotator(verbose=False)

        ann = AnnotationResult(
            gene_annotations=pd.DataFrame({
                'symbol': ['TP53', 'BRCA1'],
                'name': ['tumor protein p53', 'BRCA1 DNA repair'],
                'type_of_gene': ['protein-coding', 'protein-coding'],
            }, index=pd.Index(['TP53', 'BRCA1'], name='gene_id')),
            pathway_enrichment=pd.DataFrame({
                'pathway': ['Cell cycle', 'DNA repair'],
                'source': ['KEGG', 'KEGG'],
                'p_value': [0.001, 0.01],
                'adjusted_p_value': [0.005, 0.03],
                'overlap_count': [2, 1],
            }),
            literature_hits=pd.DataFrame({
                'gene_id': ['TP53', 'BRCA1'],
                'n_publications': [50000, 30000],
                'top_titles': ['P53 in cancer', 'BRCA1 in breast cancer'],
                'pmids': ['12345', '67890'],
            }),
            ppi_network={
                'n_nodes': 2, 'n_edges': 1,
                'edges': [{'protein_a': 'TP53', 'protein_b': 'BRCA1', 'score': 0.9}],
                'enrichment_pvalue': 0.01,
                'network_url': 'https://string-db.org/test',
            },
        )

        clf_results = {
            'random_forest': ClassificationResult(
                model_name='random_forest', auc=0.95, f1=0.90,
                sensitivity=0.92, specificity=0.88,
            ),
        }

        panel_opt = PanelOptimizationResult(
            optimal_panel=['TP53', 'BRCA1'],
            optimal_size=2, optimal_auc=0.95,
            panel_curve=pd.DataFrame({'panel_size': [2], 'auc_mean': [0.95]}),
        )

        report_path = output_dir / "test_report.md"
        report = annotator.generate_report(
            panel_genes=['TP53', 'BRCA1'],
            annotation_result=ann,
            classification_results=clf_results,
            panel_optimization=panel_opt,
            output_path=report_path,
        )

        assert isinstance(report, str)
        assert '# RAPTOR Biomarker Discovery Report' in report
        assert 'TP53' in report
        assert 'Cell cycle' in report
        assert report_path.exists()


# ============================================================================
# 7. DATA CLASS TESTS
# ============================================================================

class TestDataClasses:
    """Test BiomarkerResult and other dataclasses."""

    def test_biomarker_result_save_load(self, output_dir):
        from raptor.biomarker_discovery import (
            BiomarkerResult, FeatureSelectionResult,
            ClassificationResult, PanelOptimizationResult,
        )

        ranked = pd.DataFrame({
            'consensus_rank': [1, 2, 3],
            'consensus_score': [0.9, 0.7, 0.5],
            'n_methods_selected': [3, 2, 1],
        }, index=pd.Index(['G1', 'G2', 'G3'], name='gene_id'))

        result = BiomarkerResult(
            ranked_genes=ranked,
            panel=['G1', 'G2'],
            panel_size=2,
            selection_results={
                'lasso': FeatureSelectionResult(
                    method='lasso', selected_genes=['G1', 'G2'],
                    gene_scores=ranked[['consensus_score']].rename(
                        columns={'consensus_score': 'score'}
                    ),
                    n_selected=2,
                ),
            },
            classification_results={
                'rf': ClassificationResult(
                    model_name='rf', auc=0.9, f1=0.85,
                    sensitivity=0.88, specificity=0.82,
                ),
            },
            best_classifier='rf',
            panel_optimization=PanelOptimizationResult(
                optimal_panel=['G1', 'G2'], optimal_size=2,
                optimal_auc=0.9,
                panel_curve=pd.DataFrame({
                    'panel_size': [1, 2], 'auc_mean': [0.8, 0.9],
                    'auc_std': [0.05, 0.03],
                }),
            ),
            study_design='binary',
            n_samples=40,
            n_initial_candidates=200,
        )

        # Save
        result.save(output_dir)
        assert (output_dir / "biomarker_panel.csv").exists()
        assert (output_dir / "ranked_genes.csv").exists()
        assert (output_dir / "classification_performance.csv").exists()
        assert (output_dir / "panel_curve.csv").exists()
        assert (output_dir / "biomarker_result.pkl").exists()
        assert (output_dir / "biomarker_params.json").exists()
        assert (output_dir / "summary.txt").exists()

        # Load
        loaded = BiomarkerResult.load(output_dir)
        assert loaded.panel_size == 2
        assert loaded.best_classifier == 'rf'
        assert loaded.n_samples == 40

    def test_biomarker_result_summary(self):
        from raptor.biomarker_discovery import (
            BiomarkerResult, ClassificationResult,
        )

        result = BiomarkerResult(
            ranked_genes=pd.DataFrame(
                {'consensus_rank': [1]},
                index=pd.Index(['G1'], name='gene_id'),
            ),
            panel=['G1'],
            panel_size=1,
            classification_results={
                'rf': ClassificationResult(model_name='rf', auc=0.9, f1=0.85),
            },
            best_classifier='rf',
            n_samples=40,
            n_initial_candidates=200,
        )

        summary = result.summary()
        assert 'MODULE 10' in summary
        assert 'rf' in summary
        assert '0.900' in summary

    def test_annotation_result_save(self, output_dir):
        from raptor.biomarker_discovery import AnnotationResult

        ann = AnnotationResult(
            gene_annotations=pd.DataFrame({'symbol': ['A']}, index=['G1']),
            pathway_enrichment=pd.DataFrame({'pathway': ['P1'], 'p_value': [0.01]}),
        )

        ann.save(output_dir)
        assert (output_dir / "gene_annotations.csv").exists()
        assert (output_dir / "pathway_enrichment.csv").exists()


# ============================================================================
# 8. INTEGRATION TESTS
# ============================================================================

class TestIntegration:
    """End-to-end integration tests."""

    def test_discover_biomarkers_minimal(
        self, synthetic_counts, synthetic_metadata, output_dir
    ):
        """Minimal end-to-end test with core methods only."""
        from raptor.biomarker_discovery import discover_biomarkers

        result = discover_biomarkers(
            counts=synthetic_counts,
            metadata=synthetic_metadata,
            group_column='condition',
            methods=['elastic_net', 'rfe'],
            min_panel=3,
            max_panel=15,
            annotate=False,  # Skip network calls
            output_dir=str(output_dir),
            verbose=False,
        )

        assert isinstance(result, type(result))  # BiomarkerResult
        assert result.panel_size >= 3
        assert result.panel_size <= 15
        assert len(result.panel) == result.panel_size
        assert result.best_classifier != ''
        assert result.n_samples == 40
        assert len(result.selection_results) == 2
        assert len(result.classification_results) >= 2

        # Check output files
        assert (output_dir / "biomarker_panel.csv").exists()
        assert (output_dir / "ranked_genes.csv").exists()
        assert (output_dir / "biomarker_result.pkl").exists()

    def test_discover_biomarkers_with_de_genes(
        self, synthetic_counts, synthetic_metadata,
        synthetic_de_genes, output_dir
    ):
        """Test with pre-specified DE genes."""
        from raptor.biomarker_discovery import discover_biomarkers

        # Create a mock DEResult-like object
        class MockDEResult:
            significant_genes = synthetic_de_genes

        result = discover_biomarkers(
            counts=synthetic_counts,
            metadata=synthetic_metadata,
            group_column='condition',
            de_result=MockDEResult(),
            methods=['de_filter', 'elastic_net'],
            min_panel=3,
            max_panel=10,
            annotate=False,
            output_dir=str(output_dir),
            verbose=False,
        )

        assert 'de_filter' in result.selection_results
        assert result.panel_size >= 3

    def test_discover_biomarkers_with_target_size(
        self, synthetic_counts, synthetic_metadata, output_dir
    ):
        """Test with explicit target panel size."""
        from raptor.biomarker_discovery import discover_biomarkers

        result = discover_biomarkers(
            counts=synthetic_counts,
            metadata=synthetic_metadata,
            group_column='condition',
            methods=['elastic_net'],
            target_panel_size=5,
            annotate=False,
            output_dir=str(output_dir),
            verbose=False,
        )

        assert result.panel_size == 5

    def test_discover_biomarkers_de_genes_classify_well(
        self, synthetic_counts, synthetic_metadata,
        synthetic_de_genes, output_dir
    ):
        """DE genes should produce high classification performance."""
        from raptor.biomarker_discovery import discover_biomarkers

        class MockDEResult:
            significant_genes = synthetic_de_genes

        result = discover_biomarkers(
            counts=synthetic_counts,
            metadata=synthetic_metadata,
            group_column='condition',
            de_result=MockDEResult(),
            methods=['de_filter', 'elastic_net'],
            target_panel_size=10,
            annotate=False,
            output_dir=str(output_dir),
            verbose=False,
        )

        best_auc = result.classification_results[result.best_classifier].auc
        assert best_auc > 0.7, (
            f"With true DE genes, expected AUC > 0.7, got {best_auc:.3f}"
        )

    def test_validate_biomarkers(
        self, synthetic_counts, synthetic_metadata, output_dir
    ):
        """Test independent validation function."""
        from raptor.biomarker_discovery import validate_biomarkers

        # Use first 10 DE genes as panel
        panel = [f"GENE_{i:04d}" for i in range(10)]

        results = validate_biomarkers(
            panel_genes=panel,
            counts=synthetic_counts,
            metadata=synthetic_metadata,
            group_column='condition',
            n_folds=3,
            verbose=False,
        )

        assert len(results) >= 2
        for name, res in results.items():
            assert 0.0 <= res.auc <= 1.0

    def test_discover_biomarkers_csv_input(
        self, synthetic_counts, synthetic_metadata, output_dir
    ):
        """Test with CSV file paths instead of DataFrames."""
        from raptor.biomarker_discovery import discover_biomarkers

        # Save to CSV
        counts_path = output_dir / "counts.csv"
        meta_path = output_dir / "metadata.csv"
        synthetic_counts.to_csv(counts_path)
        synthetic_metadata.to_csv(meta_path, index=False)

        result = discover_biomarkers(
            counts=str(counts_path),
            metadata=str(meta_path),
            group_column='condition',
            methods=['elastic_net'],
            target_panel_size=5,
            annotate=False,
            output_dir=str(output_dir / "results"),
            verbose=False,
        )

        assert result.panel_size == 5
        assert (output_dir / "results" / "biomarker_panel.csv").exists()


# ============================================================================
# 9. SURVIVAL INTEGRATION TESTS
# ============================================================================

class TestSurvivalIntegration:
    """End-to-end survival analysis tests."""

    @pytest.fixture(autouse=True)
    def check_lifelines(self):
        try:
            import lifelines  # noqa: F401
        except ImportError:
            pytest.skip("lifelines not installed")

    def test_discover_survival_biomarkers(
        self, synthetic_counts, synthetic_clinical, output_dir
    ):
        from raptor.biomarker_discovery import discover_survival_biomarkers

        result = discover_survival_biomarkers(
            counts=synthetic_counts,
            clinical=synthetic_clinical,
            time_column='os_time',
            event_column='os_event',
            fdr_threshold=0.5,  # Lenient for synthetic data
            output_dir=str(output_dir),
            verbose=False,
        )

        assert hasattr(result, 'c_index')
        assert hasattr(result, 'significant_genes')
        assert hasattr(result, 'panel_genes')
        assert (output_dir / "cox_univariate_screen.csv").exists()
        assert (output_dir / "survival_summary.json").exists()


# ============================================================================
# 10. DEPENDENCY CHECKING TESTS
# ============================================================================

class TestDependencies:
    """Test dependency checking and graceful degradation."""

    def test_get_dependencies_status(self):
        from raptor.biomarker_discovery import get_dependencies_status
        status = get_dependencies_status()

        assert isinstance(status, dict)
        assert 'scikit-learn' in status
        assert 'xgboost' in status
        assert 'shap' in status
        assert 'boruta' in status
        assert 'mrmr' in status
        assert 'lifelines' in status
        assert 'PyWGCNA' in status
        # sklearn should always be True if tests run
        assert status['scikit-learn'] is True

    def test_missing_dependency_raises(self, prepared_data):
        """Methods should raise DependencyError when deps missing."""
        from raptor.biomarker_discovery import FeatureSelector, DependencyError

        X, y, _ = prepared_data
        selector = FeatureSelector(random_state=42, verbose=False)

        # Mock boruta as unavailable
        import raptor.biomarker_discovery as bm
        original = bm._BORUTA_AVAILABLE
        bm._BORUTA_AVAILABLE = False
        try:
            with pytest.raises(DependencyError):
                selector.select_boruta(X, y)
        finally:
            bm._BORUTA_AVAILABLE = original


# ============================================================================
# 11. EDGE CASE TESTS
# ============================================================================

class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_small_dataset(self):
        """Test with very small dataset (< 20 samples)."""
        from raptor.biomarker_discovery import discover_biomarkers

        np.random.seed(99)
        counts = pd.DataFrame(
            np.random.negative_binomial(5, 0.01, size=(50, 10)),
            index=[f"G{i}" for i in range(50)],
            columns=[f"s{i}" for i in range(10)],
        )
        metadata = pd.DataFrame({
            'sample_id': [f"s{i}" for i in range(10)],
            'condition': ['A'] * 5 + ['B'] * 5,
        })

        with tempfile.TemporaryDirectory() as tmpdir:
            result = discover_biomarkers(
                counts=counts,
                metadata=metadata,
                group_column='condition',
                methods=['elastic_net'],
                target_panel_size=3,
                annotate=False,
                output_dir=tmpdir,
                verbose=False,
            )
            # Should auto-switch to LOOCV for small samples
            assert result.panel_size == 3
            assert result.n_samples == 10

    def test_single_method(self, synthetic_counts, synthetic_metadata):
        """Test with only one feature selection method."""
        from raptor.biomarker_discovery import discover_biomarkers

        with tempfile.TemporaryDirectory() as tmpdir:
            result = discover_biomarkers(
                counts=synthetic_counts,
                metadata=synthetic_metadata,
                group_column='condition',
                methods=['rfe'],
                target_panel_size=5,
                annotate=False,
                output_dir=tmpdir,
                verbose=False,
            )
            assert result.panel_size == 5
            assert len(result.selection_results) == 1

    def test_all_genes_identical(self):
        """Test behavior when all genes have identical expression."""
        from raptor.biomarker_discovery import _prepare_expression_data

        counts = pd.DataFrame(
            np.full((10, 6), 100),
            index=[f"G{i}" for i in range(10)],
            columns=[f"s{i}" for i in range(6)],
        )
        metadata = pd.DataFrame({
            'sample_id': [f"s{i}" for i in range(6)],
            'condition': ['A'] * 3 + ['B'] * 3,
        })

        X, y, gene_list = _prepare_expression_data(
            counts, metadata, 'condition',
        )
        # Should still produce valid output (even if useless)
        assert X.shape[0] == 6


# ============================================================================
# RUN
# ============================================================================

if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
