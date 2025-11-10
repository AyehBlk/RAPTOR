#!/usr/bin/env python3

"""
02_profile_data.py
Profile RNA-seq count data to characterize dataset properties

Usage: python 02_profile_data.py <count_matrix.csv> [--output profile.json]

Author: Ayeh Bolouki
Organization: RAPTOR Project
License: MIT
"""

import pandas as pd
import numpy as np
import argparse
import json
import sys
from pathlib import Path
from scipy import stats

class RNAseqDataProfiler:
    """Profile RNA-seq count data characteristics"""
    
    def __init__(self, count_matrix):
        """
        Initialize profiler with count matrix
        
        Parameters
        ----------
        count_matrix : pd.DataFrame
            Count matrix with genes as rows and samples as columns
        """
        self.counts = count_matrix
        self.profile = {}
        
    def calculate_basic_stats(self):
        """Calculate basic dataset statistics"""
        print("Calculating basic statistics...")
        
        self.profile['n_genes'] = len(self.counts)
        self.profile['n_samples'] = len(self.counts.columns)
        
        # Library sizes
        library_sizes = self.counts.sum(axis=0)
        self.profile['mean_library_size'] = float(library_sizes.mean())
        self.profile['median_library_size'] = float(library_sizes.median())
        self.profile['library_size_cv'] = float(library_sizes.std() / library_sizes.mean())
        
        # Gene detection
        genes_detected = (self.counts > 0).sum(axis=1)
        self.profile['mean_genes_detected'] = float(genes_detected.mean())
        
        # Zero inflation
        total_values = self.counts.size
        zero_values = (self.counts == 0).sum().sum()
        self.profile['zero_percentage'] = float(100 * zero_values / total_values)
        
        return self.profile
    
    def calculate_expression_characteristics(self):
        """Calculate expression-level characteristics"""
        print("Analyzing expression characteristics...")
        
        # Mean expression per gene
        gene_means = self.counts.mean(axis=1)
        self.profile['mean_expression'] = float(gene_means.mean())
        self.profile['median_expression'] = float(gene_means.median())
        
        # Expression variance
        gene_vars = self.counts.var(axis=1)
        self.profile['mean_variance'] = float(gene_vars.mean())
        
        # Low vs high expression genes
        low_expr = (gene_means < 10).sum()
        high_expr = (gene_means > 100).sum()
        self.profile['low_expression_genes'] = int(low_expr)
        self.profile['high_expression_genes'] = int(high_expr)
        self.profile['pct_low_expression'] = float(100 * low_expr / len(gene_means))
        self.profile['pct_high_expression'] = float(100 * high_expr / len(gene_means))
        
        # Dynamic range
        nonzero_means = gene_means[gene_means > 0]
        if len(nonzero_means) > 0:
            self.profile['expression_range_log10'] = float(
                np.log10(nonzero_means.max()) - np.log10(nonzero_means.min())
            )
        
        return self.profile
    
    def calculate_dispersion(self):
        """Calculate biological coefficient of variation"""
        print("Calculating dispersion...")
        
        gene_means = self.counts.mean(axis=1)
        gene_vars = self.counts.var(axis=1)
        
        # Biological coefficient of variation (BCV)
        # BCV^2 = (Var - Mean) / Mean^2 for genes with mean > 0
        mask = gene_means > 0
        bcv_squared = (gene_vars[mask] - gene_means[mask]) / (gene_means[mask] ** 2)
        bcv_squared = bcv_squared[bcv_squared > 0]  # Remove negative values
        
        if len(bcv_squared) > 0:
            self.profile['median_bcv'] = float(np.sqrt(bcv_squared.median()))
            self.profile['mean_bcv'] = float(np.sqrt(bcv_squared.mean()))
        
        return self.profile
    
    def calculate_sample_correlation(self):
        """Calculate inter-sample correlations"""
        print("Calculating sample correlations...")
        
        # Log-transform for correlation
        log_counts = np.log2(self.counts + 1)
        
        # Pairwise correlations
        corr_matrix = log_counts.corr()
        
        # Upper triangle (excluding diagonal)
        upper_tri = corr_matrix.where(
            np.triu(np.ones(corr_matrix.shape), k=1).astype(bool)
        )
        
        correlations = upper_tri.stack().values
        self.profile['mean_sample_correlation'] = float(correlations.mean())
        self.profile['median_sample_correlation'] = float(np.median(correlations))
        self.profile['min_sample_correlation'] = float(correlations.min())
        
        return self.profile
    
    def calculate_complexity(self):
        """Calculate library complexity metrics"""
        print("Calculating library complexity...")
        
        # Genes detected per sample
        genes_per_sample = (self.counts > 0).sum(axis=0)
        self.profile['mean_genes_per_sample'] = float(genes_per_sample.mean())
        self.profile['min_genes_per_sample'] = int(genes_per_sample.min())
        
        # Duplication rate estimate (genes with very high counts)
        library_sizes = self.counts.sum(axis=0)
        top_gene_pct = []
        
        for col in self.counts.columns:
            sorted_counts = self.counts[col].sort_values(ascending=False)
            top_20_sum = sorted_counts.head(20).sum()
            top_gene_pct.append(100 * top_20_sum / sorted_counts.sum())
        
        self.profile['mean_top20_percentage'] = float(np.mean(top_gene_pct))
        
        return self.profile
    
    def detect_batch_effects(self):
        """Simple batch effect detection"""
        print("Checking for potential batch effects...")
        
        # PCA to check for outliers
        from sklearn.decomposition import PCA
        
        log_counts = np.log2(self.counts + 1).T
        
        try:
            pca = PCA(n_components=min(3, len(self.counts.columns)))
            pca_result = pca.fit_transform(log_counts)
            
            self.profile['pca_var_explained_pc1'] = float(pca.explained_variance_ratio_[0])
            self.profile['pca_var_explained_pc2'] = float(pca.explained_variance_ratio_[1]) if len(pca.explained_variance_ratio_) > 1 else 0
            
            # Check for outliers (samples > 3 SD from mean on PC1)
            pc1_scores = pca_result[:, 0]
            outliers = np.abs(stats.zscore(pc1_scores)) > 3
            self.profile['potential_outliers'] = int(outliers.sum())
            
        except Exception as e:
            print(f"  Warning: PCA analysis failed: {e}")
            self.profile['pca_var_explained_pc1'] = None
        
        return self.profile
    
    def generate_full_profile(self):
        """Generate complete data profile"""
        print("\n=== RNA-seq Data Profiling ===")
        print(f"Samples: {len(self.counts.columns)}")
        print(f"Genes: {len(self.counts)}")
        print()
        
        self.calculate_basic_stats()
        self.calculate_expression_characteristics()
        self.calculate_dispersion()
        self.calculate_sample_correlation()
        self.calculate_complexity()
        self.detect_batch_effects()
        
        return self.profile
    
    def recommend_pipeline(self):
        """Recommend suitable pipelines based on profile"""
        print("\n=== Pipeline Recommendations ===")
        
        recommendations = []
        
        # Based on sample size
        n_samples = self.profile['n_samples']
        if n_samples < 2:
            recommendations.append("Pipeline 6 (NOISeq) - designed for no replicates")
        elif n_samples <= 3:
            recommendations.append("Pipeline 7 (EBSeq) - optimized for small sample sizes")
        elif n_samples >= 6:
            recommendations.append("Pipeline 1 (STAR-RSEM-DESeq2) - gold standard for adequate samples")
            recommendations.append("Pipeline 5 (limma-voom) - excellent for complex designs")
        
        # Based on library size
        if self.profile['mean_library_size'] < 10e6:
            recommendations.append("Pipeline 3 (Salmon) or Pipeline 4 (Kallisto) - efficient for lower coverage")
        
        # Based on zero inflation
        if self.profile['zero_percentage'] > 70:
            recommendations.append("Consider filtering lowly expressed genes")
        
        # Based on correlation (batch effects)
        if self.profile['min_sample_correlation'] < 0.7:
            recommendations.append("Pipeline 5 (limma-voom) - can handle batch effects")
            recommendations.append("Consider batch correction")
        
        # General speed recommendation
        if n_samples > 20:
            recommendations.append("Pipeline 3 (Salmon) or Pipeline 4 (Kallisto) - faster for large datasets")
        
        self.profile['recommendations'] = recommendations
        
        for i, rec in enumerate(recommendations, 1):
            print(f"{i}. {rec}")
        
        return recommendations


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description='Profile RNA-seq count data for RAPTOR pipeline selection'
    )
    parser.add_argument('count_matrix', 
                       help='Path to count matrix CSV file (genes x samples)')
    parser.add_argument('-o', '--output', 
                       default='data_profile.json',
                       help='Output JSON file for profile (default: data_profile.json)')
    parser.add_argument('-v', '--verbose', 
                       action='store_true',
                       help='Verbose output')
    
    args = parser.parse_args()
    
    # Load count matrix
    print(f"Loading count matrix from: {args.count_matrix}")
    try:
        counts = pd.read_csv(args.count_matrix, index_col=0)
    except Exception as e:
        print(f"Error loading count matrix: {e}")
        sys.exit(1)
    
    # Create profiler
    profiler = RNAseqDataProfiler(counts)
    
    # Generate profile
    profile = profiler.generate_full_profile()
    
    # Get recommendations
    profiler.recommend_pipeline()
    
    # Save profile
    output_path = Path(args.output)
    with open(output_path, 'w') as f:
        json.dump(profile, f, indent=2)
    
    print(f"\n✓ Profile saved to: {output_path}")
    
    # Print summary
    print("\n=== Profile Summary ===")
    print(f"Samples: {profile['n_samples']}")
    print(f"Genes: {profile['n_genes']}")
    print(f"Mean library size: {profile['mean_library_size']:,.0f}")
    print(f"Library size CV: {profile['library_size_cv']:.2f}")
    print(f"Zero percentage: {profile['zero_percentage']:.1f}%")
    print(f"Mean sample correlation: {profile['mean_sample_correlation']:.3f}")
    
    if 'median_bcv' in profile:
        print(f"Median BCV: {profile['median_bcv']:.3f}")
    
    if profile.get('potential_outliers', 0) > 0:
        print(f"\n⚠ Warning: {profile['potential_outliers']} potential outlier sample(s) detected")
    
    print("\n✓ Data profiling complete!")


if __name__ == '__main__':
    main()
