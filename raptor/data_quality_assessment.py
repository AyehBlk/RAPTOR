#!/usr/bin/env python3

"""
Advanced Data Quality Assessment Module

Provides comprehensive quality metrics including:
- Batch effect detection
- Overall quality scoring
- Sample clustering analysis
- Outlier detection
- Technical variation assessment
- Biological signal strength

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

import logging
logger = logging.getLogger(__name__)


class DataQualityAssessor:
    """
    Comprehensive data quality assessment for RNA-seq data.
    
    Performs:
    - Batch effect detection
    - Quality scoring
    - Outlier detection
    - Sample clustering
    - Technical variation assessment
    
    Parameters
    ----------
    counts : pd.DataFrame
        Count matrix (genes x samples)
    metadata : pd.DataFrame, optional
        Sample metadata with batch information
    
    Examples
    --------
    >>> assessor = DataQualityAssessor(counts, metadata)
    >>> quality_report = assessor.assess_quality()
    >>> print(f"Overall Quality Score: {quality_report['overall_score']:.2f}")
    """
    
    def __init__(self, counts, metadata=None):
        """Initialize quality assessor."""
        self.counts = counts
        self.metadata = metadata
        self.n_genes, self.n_samples = counts.shape
        
        # Results storage
        self.results = {}
        
        logger.info(f"Initialized quality assessor for {self.n_samples} samples, {self.n_genes} genes")
    
    def assess_quality(self):
        """
        Perform comprehensive quality assessment.
        
        Returns
        -------
        dict
            Complete quality report with all metrics
        """
        logger.info("Starting comprehensive quality assessment...")
        
        # Run all assessments
        self._assess_library_quality()
        self._assess_gene_detection()
        self._detect_outliers()
        self._assess_variance_structure()
        self._detect_batch_effects()
        self._assess_biological_signal()
        self._calculate_overall_score()
        
        # Compile report
        report = self._compile_report()
        
        logger.info("Quality assessment complete")
        return report
    
    def _assess_library_quality(self):
        """Assess library size and distribution quality."""
        logger.info("Assessing library quality...")
        
        lib_sizes = self.counts.sum(axis=0)
        
        # Library size metrics
        mean_size = lib_sizes.mean()
        cv = lib_sizes.std() / mean_size
        min_size = lib_sizes.min()
        max_size = lib_sizes.max()
        
        # Quality flags
        flags = []
        if cv > 0.5:
            flags.append("High library size variation")
        if min_size < mean_size * 0.3:
            flags.append("Some libraries very small")
        if max_size > mean_size * 3:
            flags.append("Some libraries very large")
        
        # Quality score (0-100)
        score = 100
        score -= min(30, cv * 50)  # Penalize high CV
        score -= min(20, (max_size / min_size - 1) * 10)  # Penalize size range
        
        self.results['library_quality'] = {
            'mean_size': float(mean_size),
            'cv': float(cv),
            'min_size': float(min_size),
            'max_size': float(max_size),
            'size_range': float(max_size / min_size),
            'score': max(0, score),
            'flags': flags,
            'status': 'good' if score > 70 else 'warning' if score > 50 else 'poor'
        }
    
    def _assess_gene_detection(self):
        """Assess gene detection and expression distribution."""
        logger.info("Assessing gene detection...")
        
        # Detection metrics
        zero_counts = (self.counts == 0).sum(axis=1)
        detection_rate = 1 - (zero_counts / self.n_samples)
        
        # Expression levels
        mean_expr = self.counts.mean(axis=1)
        median_expr = self.counts.median(axis=1)
        
        # Categories
        n_high = (mean_expr > mean_expr.quantile(0.75)).sum()
        n_medium = ((mean_expr > mean_expr.quantile(0.25)) & 
                    (mean_expr <= mean_expr.quantile(0.75))).sum()
        n_low = (mean_expr <= mean_expr.quantile(0.25)).sum()
        
        # Quality flags
        flags = []
        zero_pct = (zero_counts.mean() / self.n_samples) * 100
        if zero_pct > 70:
            flags.append("Very high zero inflation")
        elif zero_pct > 50:
            flags.append("High zero inflation")
        
        # Quality score
        score = 100
        score -= min(40, zero_pct * 0.5)  # Penalize zeros
        if n_high < self.n_genes * 0.05:
            score -= 10  # Few highly expressed genes
        
        self.results['gene_detection'] = {
            'mean_detection_rate': float(detection_rate.mean()),
            'zero_inflation_pct': float(zero_pct),
            'n_highly_expressed': int(n_high),
            'n_medium_expressed': int(n_medium),
            'n_low_expressed': int(n_low),
            'score': max(0, score),
            'flags': flags,
            'status': 'good' if score > 70 else 'warning' if score > 50 else 'poor'
        }
    
    def _detect_outliers(self):
        """Detect outlier samples using multiple methods."""
        logger.info("Detecting outlier samples...")
        
        # Normalize for outlier detection
        normalized = np.log2(self.counts + 1)
        
        # Method 1: PCA-based outlier detection
        pca = PCA(n_components=min(5, self.n_samples - 1))
        pca_scores = pca.fit_transform(normalized.T)
        
        # Mahalanobis distance
        mean = pca_scores.mean(axis=0)
        cov = np.cov(pca_scores.T)
        inv_cov = np.linalg.pinv(cov)
        
        mahal_dist = []
        for i in range(len(pca_scores)):
            diff = pca_scores[i] - mean
            dist = np.sqrt(diff @ inv_cov @ diff.T)
            mahal_dist.append(dist)
        
        mahal_dist = np.array(mahal_dist)
        
        # Outlier threshold (3 standard deviations)
        threshold = mahal_dist.mean() + 3 * mahal_dist.std()
        outliers_pca = mahal_dist > threshold
        
        # Method 2: Hierarchical clustering
        distances = pdist(normalized.T, metric='correlation')
        linkage_matrix = linkage(distances, method='average')
        
        # Method 3: Library size outliers
        lib_sizes = self.counts.sum(axis=0)
        z_scores = np.abs(stats.zscore(lib_sizes))
        outliers_libsize = z_scores > 3
        
        # Combine methods
        outlier_scores = (outliers_pca.astype(int) + 
                         outliers_libsize.astype(int))
        
        outlier_samples = self.counts.columns[outlier_scores >= 2].tolist()
        
        # Quality score
        n_outliers = len(outlier_samples)
        outlier_pct = (n_outliers / self.n_samples) * 100
        
        score = 100 - min(50, outlier_pct * 5)
        
        flags = []
        if n_outliers > 0:
            flags.append(f"{n_outliers} potential outlier sample(s) detected")
        
        self.results['outlier_detection'] = {
            'n_outliers': n_outliers,
            'outlier_percentage': float(outlier_pct),
            'outlier_samples': outlier_samples,
            'mahalanobis_distances': mahal_dist.tolist(),
            'score': max(0, score),
            'flags': flags,
            'status': 'good' if n_outliers == 0 else 'warning' if n_outliers <= 2 else 'poor'
        }
    
    def _assess_variance_structure(self):
        """Assess variance structure and components."""
        logger.info("Assessing variance structure...")
        
        # Normalize
        normalized = np.log2(self.counts + 1)
        
        # PCA analysis
        scaler = StandardScaler()
        scaled = scaler.fit_transform(normalized.T)
        
        pca = PCA()
        pca.fit(scaled)
        
        variance_explained = pca.explained_variance_ratio_
        
        # First PC dominance
        pc1_var = variance_explained[0] * 100
        pc2_var = variance_explained[1] * 100 if len(variance_explained) > 1 else 0
        
        # Cumulative variance
        cum_var = np.cumsum(variance_explained) * 100
        n_components_90 = np.where(cum_var >= 90)[0][0] + 1 if any(cum_var >= 90) else len(variance_explained)
        
        # Quality assessment
        flags = []
        if pc1_var > 50:
            flags.append("First PC explains >50% variance (potential batch effect)")
        if n_components_90 < 3:
            flags.append("Very few components needed for 90% variance")
        
        # Score
        score = 100
        if pc1_var > 50:
            score -= 20
        if pc1_var > 70:
            score -= 20
        if n_components_90 < 2:
            score -= 15
        
        self.results['variance_structure'] = {
            'pc1_variance': float(pc1_var),
            'pc2_variance': float(pc2_var),
            'variance_explained_top5': variance_explained[:5].tolist(),
            'n_components_90pct': int(n_components_90),
            'score': max(0, score),
            'flags': flags,
            'status': 'good' if score > 70 else 'warning' if score > 50 else 'poor'
        }
    
    def _detect_batch_effects(self):
        """Detect potential batch effects."""
        logger.info("Detecting batch effects...")
        
        # Normalize
        normalized = np.log2(self.counts + 1)
        
        # PCA for visualization
        pca = PCA(n_components=min(3, self.n_samples - 1))
        pca_coords = pca.fit_transform(normalized.T)
        
        batch_detected = False
        batch_strength = 0
        flags = []
        batch_variable = None
        
        if self.metadata is not None and len(self.metadata) > 0:
            # Check each categorical variable for batch effects
            categorical_cols = self.metadata.select_dtypes(include=['object', 'category']).columns
            
            batch_scores = {}
            
            for col in categorical_cols:
                if col in self.metadata.columns:
                    groups = self.metadata[col].values
                    unique_groups = np.unique(groups)
                    
                    if len(unique_groups) > 1 and len(unique_groups) < self.n_samples:
                        # Perform PERMANOVA-like test on PC1
                        try:
                            # Calculate between-group vs within-group variance
                            between_var = 0
                            within_var = 0
                            
                            overall_mean = pca_coords[:, 0].mean()
                            
                            for group in unique_groups:
                                group_mask = groups == group
                                group_data = pca_coords[group_mask, 0]
                                
                                if len(group_data) > 0:
                                    group_mean = group_data.mean()
                                    between_var += len(group_data) * (group_mean - overall_mean) ** 2
                                    within_var += ((group_data - group_mean) ** 2).sum()
                            
                            if within_var > 0:
                                f_stat = (between_var / (len(unique_groups) - 1)) / (within_var / (self.n_samples - len(unique_groups)))
                                batch_scores[col] = f_stat
                        except:
                            pass
            
            if batch_scores:
                # Find strongest batch effect
                batch_variable = max(batch_scores.items(), key=lambda x: x[1])[0]
                batch_strength = batch_scores[batch_variable]
                
                if batch_strength > 5:  # Arbitrary threshold
                    batch_detected = True
                    flags.append(f"Strong batch effect detected in '{batch_variable}'")
                elif batch_strength > 2:
                    flags.append(f"Moderate batch effect possible in '{batch_variable}'")
        else:
            flags.append("No metadata provided - cannot check for batch effects")
        
        # Clustering-based detection (without metadata)
        # Check if samples cluster into distinct groups
        from scipy.cluster.hierarchy import fcluster
        distances = pdist(normalized.T, metric='correlation')
        linkage_matrix = linkage(distances, method='ward')
        
        # Try different numbers of clusters
        best_silhouette = -1
        for n_clusters in range(2, min(6, self.n_samples)):
            clusters = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
            
            # Calculate silhouette-like score
            cluster_separation = 0
            for i in range(n_clusters):
                cluster_mask = clusters == (i + 1)
                if cluster_mask.sum() > 1:
                    cluster_data = pca_coords[cluster_mask, 0]
                    cluster_mean = cluster_data.mean()
                    cluster_std = cluster_data.std()
                    if cluster_std > 0:
                        separation = abs(cluster_mean) / cluster_std
                        cluster_separation += separation
            
            avg_separation = cluster_separation / n_clusters
            if avg_separation > best_silhouette:
                best_silhouette = avg_separation
        
        if best_silhouette > 1.5 and not batch_detected:
            flags.append("Samples form distinct clusters (possible hidden batch effect)")
        
        # Score
        score = 100
        if batch_detected:
            score -= min(40, batch_strength * 5)
        if best_silhouette > 2:
            score -= 15
        
        self.results['batch_effects'] = {
            'batch_detected': batch_detected,
            'batch_variable': batch_variable,
            'batch_strength': float(batch_strength) if batch_strength > 0 else None,
            'cluster_separation': float(best_silhouette),
            'pca_coordinates': pca_coords.tolist(),
            'score': max(0, score),
            'flags': flags,
            'status': 'good' if score > 80 else 'warning' if score > 60 else 'poor',
            'recommendation': self._get_batch_recommendation(batch_detected, batch_strength)
        }
    
    def _get_batch_recommendation(self, detected, strength):
        """Get recommendation for batch effect handling."""
        if not detected:
            return "No batch correction needed"
        
        if strength > 10:
            return "Strong batch effect - Use ComBat or limma removeBatchEffect before analysis"
        elif strength > 5:
            return "Moderate batch effect - Consider including batch as covariate in DE analysis"
        else:
            return "Weak batch effect - Monitor in exploratory analysis"
    
    def _assess_biological_signal(self):
        """Assess strength of biological signal."""
        logger.info("Assessing biological signal strength...")
        
        # Calculate coefficient of variation for each gene
        means = self.counts.mean(axis=1)
        stds = self.counts.std(axis=1)
        
        # Avoid division by zero
        cv = np.zeros(len(means))
        nonzero = means > 0
        cv[nonzero] = stds[nonzero] / means[nonzero]
        
        # Biological vs technical variation
        # Genes with high variance relative to mean likely biological
        median_cv = np.median(cv[cv > 0])
        
        # Count highly variable genes
        high_var_threshold = np.percentile(cv, 75)
        n_high_var = (cv > high_var_threshold).sum()
        
        # Signal-to-noise estimate
        # High mean, high variance = good signal
        high_mean = means > means.quantile(0.5)
        high_var = cv > median_cv
        signal_genes = high_mean & high_var
        n_signal = signal_genes.sum()
        
        signal_strength = (n_signal / self.n_genes) * 100
        
        flags = []
        if signal_strength < 10:
            flags.append("Low biological signal detected")
        elif signal_strength > 30:
            flags.append("Strong biological signal detected")
        
        # Score
        score = min(100, signal_strength * 2)
        
        self.results['biological_signal'] = {
            'median_cv': float(median_cv),
            'n_high_variance_genes': int(n_high_var),
            'n_signal_genes': int(n_signal),
            'signal_strength_pct': float(signal_strength),
            'score': score,
            'flags': flags,
            'status': 'good' if score > 60 else 'warning' if score > 40 else 'poor'
        }
    
    def _calculate_overall_score(self):
        """Calculate overall quality score."""
        logger.info("Calculating overall quality score...")
        
        # Weighted average of component scores
        weights = {
            'library_quality': 0.20,
            'gene_detection': 0.20,
            'outlier_detection': 0.15,
            'variance_structure': 0.15,
            'batch_effects': 0.15,
            'biological_signal': 0.15
        }
        
        overall = 0
        for component, weight in weights.items():
            if component in self.results:
                overall += self.results[component]['score'] * weight
        
        # Overall status
        if overall >= 80:
            status = 'excellent'
            recommendation = "High quality data - proceed with analysis"
        elif overall >= 70:
            status = 'good'
            recommendation = "Good quality data - suitable for analysis"
        elif overall >= 60:
            status = 'acceptable'
            recommendation = "Acceptable quality - review flagged issues"
        elif overall >= 50:
            status = 'poor'
            recommendation = "Poor quality - address issues before analysis"
        else:
            status = 'fail'
            recommendation = "Failed quality check - major issues need resolution"
        
        self.results['overall'] = {
            'score': overall,
            'status': status,
            'recommendation': recommendation
        }
    
    def _compile_report(self):
        """Compile comprehensive quality report."""
        report = {
            'dataset_info': {
                'n_samples': self.n_samples,
                'n_genes': self.n_genes,
                'has_metadata': self.metadata is not None
            },
            'overall': self.results['overall'],
            'components': {
                'library_quality': self.results['library_quality'],
                'gene_detection': self.results['gene_detection'],
                'outlier_detection': self.results['outlier_detection'],
                'variance_structure': self.results['variance_structure'],
                'batch_effects': self.results['batch_effects'],
                'biological_signal': self.results['biological_signal']
            },
            'summary': self._generate_summary()
        }
        
        return report
    
    def _generate_summary(self):
        """Generate human-readable summary."""
        lines = []
        
        # Overall
        overall = self.results['overall']
        lines.append(f"Overall Quality: {overall['status'].upper()} ({overall['score']:.1f}/100)")
        lines.append(f"Recommendation: {overall['recommendation']}")
        lines.append("")
        
        # Component summaries
        components = [
            ('Library Quality', 'library_quality'),
            ('Gene Detection', 'gene_detection'),
            ('Outlier Detection', 'outlier_detection'),
            ('Variance Structure', 'variance_structure'),
            ('Batch Effects', 'batch_effects'),
            ('Biological Signal', 'biological_signal')
        ]
        
        for name, key in components:
            comp = self.results[key]
            status_icon = '✓' if comp['status'] == 'good' else '⚠' if comp['status'] == 'warning' else '✗'
            lines.append(f"{status_icon} {name}: {comp['status']} ({comp['score']:.1f}/100)")
            
            if comp['flags']:
                for flag in comp['flags']:
                    lines.append(f"  • {flag}")
        
        return '\n'.join(lines)
    
    def plot_quality_report(self, output_file=None):
        """
        Generate comprehensive quality visualization.
        
        Parameters
        ----------
        output_file : str, optional
            Path to save figure
        """
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        fig = plt.figure(figsize=(16, 12))
        gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
        
        # 1. Overall score gauge
        ax1 = fig.add_subplot(gs[0, 0])
        self._plot_score_gauge(ax1, self.results['overall']['score'], 'Overall Quality')
        
        # 2. Component scores
        ax2 = fig.add_subplot(gs[0, 1:])
        self._plot_component_scores(ax2)
        
        # 3. Library size distribution
        ax3 = fig.add_subplot(gs[1, 0])
        lib_sizes = self.counts.sum(axis=0)
        ax3.hist(lib_sizes / 1e6, bins=20, color='steelblue', edgecolor='black')
        ax3.set_xlabel('Library Size (millions)', fontsize=10)
        ax3.set_ylabel('Frequency', fontsize=10)
        ax3.set_title('Library Size Distribution', fontsize=12, fontweight='bold')
        ax3.axvline(lib_sizes.mean() / 1e6, color='red', linestyle='--', label='Mean')
        ax3.legend()
        
        # 4. PCA plot
        ax4 = fig.add_subplot(gs[1, 1])
        pca_coords = np.array(self.results['batch_effects']['pca_coordinates'])
        
        if self.results['outlier_detection']['outlier_samples']:
            # Color outliers differently
            outliers = self.results['outlier_detection']['outlier_samples']
            outlier_mask = np.array([s in outliers for s in self.counts.columns])
            
            ax4.scatter(pca_coords[~outlier_mask, 0], pca_coords[~outlier_mask, 1],
                       c='steelblue', s=50, alpha=0.6, label='Normal')
            ax4.scatter(pca_coords[outlier_mask, 0], pca_coords[outlier_mask, 1],
                       c='red', s=50, alpha=0.8, label='Outlier')
            ax4.legend()
        else:
            ax4.scatter(pca_coords[:, 0], pca_coords[:, 1],
                       c='steelblue', s=50, alpha=0.6)
        
        ax4.set_xlabel('PC1', fontsize=10)
        ax4.set_ylabel('PC2', fontsize=10)
        ax4.set_title('PCA Plot', fontsize=12, fontweight='bold')
        
        # 5. Variance explained
        ax5 = fig.add_subplot(gs[1, 2])
        var_exp = self.results['variance_structure']['variance_explained_top5']
        ax5.bar(range(1, len(var_exp) + 1), np.array(var_exp) * 100, color='coral')
        ax5.set_xlabel('Principal Component', fontsize=10)
        ax5.set_ylabel('Variance Explained (%)', fontsize=10)
        ax5.set_title('PCA Variance', fontsize=12, fontweight='bold')
        
        # 6. Gene detection
        ax6 = fig.add_subplot(gs[2, 0])
        detection = self.results['gene_detection']
        categories = ['High', 'Medium', 'Low']
        values = [detection['n_highly_expressed'], 
                 detection['n_medium_expressed'],
                 detection['n_low_expressed']]
        colors = ['darkgreen', 'orange', 'lightcoral']
        ax6.pie(values, labels=categories, colors=colors, autopct='%1.1f%%', startangle=90)
        ax6.set_title('Gene Expression Levels', fontsize=12, fontweight='bold')
        
        # 7. Quality status summary
        ax7 = fig.add_subplot(gs[2, 1:])
        ax7.axis('off')
        
        summary_text = self._generate_summary()
        ax7.text(0.05, 0.95, summary_text, transform=ax7.transAxes,
                fontsize=9, verticalalignment='top', family='monospace',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
        
        plt.suptitle('🦖 RAPTOR Data Quality Assessment Report',
                    fontsize=16, fontweight='bold', y=0.98)
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            logger.info(f"Quality report saved to {output_file}")
        
        return fig
    
    def _plot_score_gauge(self, ax, score, title):
        """Plot score as gauge chart."""
        # Color based on score
        if score >= 80:
            color = 'green'
        elif score >= 60:
            color = 'orange'
        else:
            color = 'red'
        
        # Draw gauge
        theta = np.linspace(0, np.pi, 100)
        r = np.ones(100)
        
        ax.plot(theta, r, 'k-', linewidth=2)
        ax.fill_between(theta, 0, r, alpha=0.1, color='gray')
        
        # Score indicator
        score_theta = (score / 100) * np.pi
        ax.plot([0, score_theta], [0, 1], color=color, linewidth=4)
        ax.plot(score_theta, 1, 'o', color=color, markersize=15)
        
        # Text
        ax.text(np.pi/2, 0.5, f'{score:.1f}', ha='center', va='center',
               fontsize=24, fontweight='bold')
        ax.text(np.pi/2, 0.2, title, ha='center', va='center',
               fontsize=12, fontweight='bold')
        
        ax.set_xlim([0, np.pi])
        ax.set_ylim([0, 1.2])
        ax.axis('off')
    
    def _plot_component_scores(self, ax):
        """Plot component scores as horizontal bar chart."""
        components = [
            'Library\nQuality',
            'Gene\nDetection',
            'Outlier\nDetection',
            'Variance\nStructure',
            'Batch\nEffects',
            'Biological\nSignal'
        ]
        
        scores = [
            self.results['library_quality']['score'],
            self.results['gene_detection']['score'],
            self.results['outlier_detection']['score'],
            self.results['variance_structure']['score'],
            self.results['batch_effects']['score'],
            self.results['biological_signal']['score']
        ]
        
        # Color by score
        colors = ['green' if s >= 80 else 'orange' if s >= 60 else 'red' for s in scores]
        
        y_pos = np.arange(len(components))
        bars = ax.barh(y_pos, scores, color=colors, alpha=0.7, edgecolor='black')
        
        # Add score labels
        for i, (bar, score) in enumerate(zip(bars, scores)):
            ax.text(score + 2, i, f'{score:.0f}', va='center', fontweight='bold')
        
        ax.set_yticks(y_pos)
        ax.set_yticklabels(components, fontsize=9)
        ax.set_xlabel('Score', fontsize=10)
        ax.set_title('Component Scores', fontsize=12, fontweight='bold')
        ax.set_xlim([0, 105])
        ax.axvline(80, color='green', linestyle='--', alpha=0.3, linewidth=1)
        ax.axvline(60, color='orange', linestyle='--', alpha=0.3, linewidth=1)
        ax.grid(axis='x', alpha=0.3)


def quick_quality_check(counts, metadata=None, plot=True, output_file='quality_report.png'):
    """
    Quick quality assessment with visualization.
    
    Parameters
    ----------
    counts : pd.DataFrame
        Count matrix
    metadata : pd.DataFrame, optional
        Sample metadata
    plot : bool
        Generate visualization
    output_file : str
        Output file for plot
    
    Returns
    -------
    dict
        Quality report
    
    Examples
    --------
    >>> report = quick_quality_check(counts, metadata)
    >>> print(report['overall']['status'])
    """
    assessor = DataQualityAssessor(counts, metadata)
    report = assessor.assess_quality()
    
    if plot:
        assessor.plot_quality_report(output_file)
    
    # Print summary
    print("\n" + "="*70)
    print("RAPTOR DATA QUALITY ASSESSMENT")
    print("="*70)
    print(report['summary'])
    print("="*70 + "\n")
    
    return report


if __name__ == '__main__':
    print("""
    RAPTOR Advanced Data Quality Assessment Module
    ==============================================
    
    Usage:
        from data_quality_assessment import DataQualityAssessor, quick_quality_check
        
        # Quick check
        report = quick_quality_check(counts, metadata)
        
        # Detailed assessment
        assessor = DataQualityAssessor(counts, metadata)
        report = assessor.assess_quality()
        assessor.plot_quality_report('quality_report.png')
    
    Features:
        ✓ Library quality assessment
        ✓ Gene detection analysis
        ✓ Outlier detection (multiple methods)
        ✓ Variance structure analysis
        ✓ Batch effect detection
        ✓ Biological signal assessment
        ✓ Overall quality scoring
        ✓ Comprehensive visualization
    """)
