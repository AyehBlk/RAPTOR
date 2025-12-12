"""
Ensemble Analysis Methods

Combines results from multiple RNA-seq pipelines to provide robust
differential expression detection through consensus approaches.

Author: Ayeh Bolouki
Email: ayeh.bolouki@unamur.be
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Union
import logging
from scipy import stats
from scipy.stats import rankdata

logger = logging.getLogger(__name__)


class EnsembleAnalyzer:
    """
    Ensemble analysis for combining pipeline results.
    
    Implements multiple ensemble methods to combine differential expression
    results from different pipelines for more robust gene detection.
    
    Parameters
    ----------
    min_pipelines : int
        Minimum number of pipelines that must agree (default: 2)
    confidence_weights : dict, optional
        Pipeline-specific confidence weights
    
    Examples
    --------
    >>> ensemble = EnsembleAnalyzer(min_pipelines=2)
    >>> results = ensemble.combine_results(
    ...     pipeline_results={
    ...         'deseq2': deseq2_df,
    ...         'edger': edger_df,
    ...         'limma': limma_df
    ...     },
    ...     method='weighted_vote'
    ... )
    >>> consensus_genes = results['consensus_genes']
    """
    
    def __init__(
        self,
        min_pipelines: int = 2,
        confidence_weights: Optional[Dict[str, float]] = None
    ):
        """Initialize ensemble analyzer."""
        self.min_pipelines = min_pipelines
        self.confidence_weights = confidence_weights or {}
        
        logger.info(f"EnsembleAnalyzer initialized: min_pipelines={min_pipelines}")
    
    def combine_results(
        self,
        pipeline_results: Dict[str, pd.DataFrame],
        method: str = 'weighted_vote',
        alpha: float = 0.05,
        lfc_threshold: float = 0.0
    ) -> Dict:
        """
        Combine results from multiple pipelines.
        
        Parameters
        ----------
        pipeline_results : dict
            Dictionary mapping pipeline names to result DataFrames.
            Each DataFrame should have columns: gene_id, log2FoldChange, pvalue, padj
        method : str
            Ensemble method: 'weighted_vote', 'rank_aggregation', 'intersection', 'union'
        alpha : float
            Significance threshold
        lfc_threshold : float
            Log2 fold-change threshold
        
        Returns
        -------
        dict
            Ensemble results with consensus genes and statistics
        """
        logger.info(f"Combining results from {len(pipeline_results)} pipelines")
        logger.info(f"Method: {method}, alpha: {alpha}, LFC threshold: {lfc_threshold}")
        
        # Validate input
        self._validate_input(pipeline_results)
        
        # Apply method
        if method == 'weighted_vote':
            result = self._weighted_vote(pipeline_results, alpha, lfc_threshold)
        elif method == 'rank_aggregation':
            result = self._rank_aggregation(pipeline_results, alpha)
        elif method == 'intersection':
            result = self._intersection(pipeline_results, alpha, lfc_threshold)
        elif method == 'union':
            result = self._union(pipeline_results, alpha, lfc_threshold)
        else:
            raise ValueError(f"Unknown method: {method}")
        
        # Calculate concordance
        result['concordance'] = self._calculate_concordance(pipeline_results, alpha)
        result['method'] = method
        result['n_pipelines'] = len(pipeline_results)
        result['pipeline_names'] = list(pipeline_results.keys())
        
        logger.info(f"Ensemble complete: {len(result['consensus_genes'])} consensus genes")
        
        return result
    
    def _validate_input(self, pipeline_results: Dict[str, pd.DataFrame]):
        """Validate pipeline results format."""
        required_cols = ['log2FoldChange', 'pvalue']
        
        for pipeline, df in pipeline_results.items():
            # Check if DataFrame
            if not isinstance(df, pd.DataFrame):
                raise ValueError(f"{pipeline}: results must be a DataFrame")
            
            # Check required columns
            missing = set(required_cols) - set(df.columns)
            if missing:
                raise ValueError(f"{pipeline}: missing columns: {missing}")
            
            # Check for gene identifiers
            if df.index.name is None and 'gene_id' not in df.columns:
                logger.warning(f"{pipeline}: no gene identifiers found, using DataFrame index")
    
    def _weighted_vote(
        self,
        pipeline_results: Dict[str, pd.DataFrame],
        alpha: float,
        lfc_threshold: float
    ) -> Dict:
        """
        Weighted voting based on pipeline confidence.
        
        Each pipeline votes for DE genes, weighted by its confidence/accuracy.
        """
        logger.info("Using weighted voting method")
        
        # Get all unique genes
        all_genes = set()
        for df in pipeline_results.values():
            all_genes.update(df.index if df.index.name else df['gene_id'])
        
        all_genes = sorted(all_genes)
        logger.info(f"Total unique genes: {len(all_genes)}")
        
        # Initialize vote matrix
        votes = pd.DataFrame(0.0, index=all_genes, columns=list(pipeline_results.keys()))
        fold_changes = pd.DataFrame(np.nan, index=all_genes, columns=list(pipeline_results.keys()))
        pvalues = pd.DataFrame(np.nan, index=all_genes, columns=list(pipeline_results.keys()))
        
        # Collect votes from each pipeline
        for pipeline, df in pipeline_results.items():
            # Get weight for this pipeline
            weight = self.confidence_weights.get(pipeline, 1.0)
            
            # Get gene index
            if df.index.name:
                gene_index = df.index
            else:
                gene_index = df['gene_id']
                df = df.set_index('gene_id')
            
            # Determine DE genes
            if 'padj' in df.columns:
                significant = (df['padj'] < alpha) & (df['log2FoldChange'].abs() > lfc_threshold)
            else:
                significant = (df['pvalue'] < alpha) & (df['log2FoldChange'].abs() > lfc_threshold)
            
            # Cast votes
            de_genes = df.index[significant]
            votes.loc[de_genes, pipeline] = weight
            
            # Store fold changes and p-values
            fold_changes.loc[df.index, pipeline] = df['log2FoldChange']
            pvalues.loc[df.index, pipeline] = df['pvalue']
        
        # Calculate total votes
        total_votes = votes.sum(axis=1)
        n_votes = (votes > 0).sum(axis=1)
        
        # Consensus: genes with votes from min_pipelines
        consensus_mask = n_votes >= self.min_pipelines
        consensus_genes = votes.index[consensus_mask].tolist()
        
        # Create consensus DataFrame
        consensus_df = pd.DataFrame({
            'gene_id': consensus_genes,
            'n_pipelines': n_votes[consensus_mask].values,
            'vote_score': total_votes[consensus_mask].values,
            'mean_log2FC': fold_changes.loc[consensus_mask].mean(axis=1).values,
            'median_log2FC': fold_changes.loc[consensus_mask].median(axis=1).values,
            'std_log2FC': fold_changes.loc[consensus_mask].std(axis=1).values,
            'min_pvalue': pvalues.loc[consensus_mask].min(axis=1).values
        })
        
        # Sort by vote score
        consensus_df = consensus_df.sort_values('vote_score', ascending=False)
        
        result = {
            'consensus_genes': consensus_df,
            'all_votes': votes,
            'fold_changes': fold_changes,
            'pvalues': pvalues
        }
        
        return result
    
    def _rank_aggregation(
        self,
        pipeline_results: Dict[str, pd.DataFrame],
        alpha: float
    ) -> Dict:
        """
        Rank aggregation using Robust Rank Aggregation (RRA).
        
        Combines gene rankings from multiple pipelines using
        order statistics.
        """
        logger.info("Using rank aggregation method (RRA)")
        
        # Get all unique genes
        all_genes = set()
        for df in pipeline_results.values():
            all_genes.update(df.index if df.index.name else df['gene_id'])
        
        all_genes = sorted(all_genes)
        n_genes = len(all_genes)
        
        # Initialize rank matrix
        ranks = pd.DataFrame(n_genes, index=all_genes, columns=list(pipeline_results.keys()))
        
        # Collect ranks from each pipeline
        for pipeline, df in pipeline_results.items():
            if df.index.name:
                gene_index = df.index
            else:
                gene_index = df['gene_id']
                df = df.set_index('gene_id')
            
            # Rank by p-value (lower is better)
            gene_ranks = rankdata(df['pvalue'], method='min')
            ranks.loc[df.index, pipeline] = gene_ranks
        
        # Calculate RRA score
        n_lists = len(pipeline_results)
        rra_scores = []
        
        for gene in all_genes:
            gene_ranks = ranks.loc[gene].values
            # RRA: beta distribution CDF
            rho = stats.beta.cdf(gene_ranks / n_genes, 1, n_lists)
            rra_score = np.min(rho) * n_lists
            rra_scores.append(rra_score)
        
        # Create result DataFrame
        rra_df = pd.DataFrame({
            'gene_id': all_genes,
            'rra_score': rra_scores,
            'mean_rank': ranks.mean(axis=1).values,
            'min_rank': ranks.min(axis=1).values,
            'n_pipelines': (ranks < n_genes).sum(axis=1).values
        })
        
        # Filter by min_pipelines
        rra_df = rra_df[rra_df['n_pipelines'] >= self.min_pipelines]
        
        # Sort by RRA score (lower is better)
        rra_df = rra_df.sort_values('rra_score')
        
        # Apply significance threshold
        consensus_genes = rra_df[rra_df['rra_score'] < alpha]
        
        result = {
            'consensus_genes': consensus_genes,
            'all_ranks': ranks,
            'rra_scores': rra_df
        }
        
        return result
    
    def _intersection(
        self,
        pipeline_results: Dict[str, pd.DataFrame],
        alpha: float,
        lfc_threshold: float
    ) -> Dict:
        """
        Intersection: genes significant in ALL pipelines.
        
        Conservative approach - only genes detected by all pipelines.
        """
        logger.info("Using intersection method (conservative)")
        
        # Get DE genes from each pipeline
        de_gene_sets = []
        
        for pipeline, df in pipeline_results.items():
            if df.index.name:
                gene_index = df.index
            else:
                gene_index = df['gene_id']
                df = df.set_index('gene_id')
            
            # Determine DE genes
            if 'padj' in df.columns:
                significant = (df['padj'] < alpha) & (df['log2FoldChange'].abs() > lfc_threshold)
            else:
                significant = (df['pvalue'] < alpha) & (df['log2FoldChange'].abs() > lfc_threshold)
            
            de_genes = set(df.index[significant])
            de_gene_sets.append(de_genes)
            logger.info(f"{pipeline}: {len(de_genes)} DE genes")
        
        # Find intersection
        consensus_genes = set.intersection(*de_gene_sets)
        logger.info(f"Intersection: {len(consensus_genes)} genes")
        
        # Create consensus DataFrame
        consensus_list = []
        for gene in consensus_genes:
            gene_data = {'gene_id': gene}
            
            for pipeline, df in pipeline_results.items():
                if df.index.name:
                    gene_df = df.loc[gene]
                else:
                    gene_df = df[df['gene_id'] == gene].iloc[0]
                
                gene_data[f'{pipeline}_log2FC'] = gene_df['log2FoldChange']
                gene_data[f'{pipeline}_pvalue'] = gene_df['pvalue']
            
            consensus_list.append(gene_data)
        
        consensus_df = pd.DataFrame(consensus_list)
        
        # Add summary statistics
        lfc_cols = [c for c in consensus_df.columns if '_log2FC' in c]
        consensus_df['mean_log2FC'] = consensus_df[lfc_cols].mean(axis=1)
        consensus_df['std_log2FC'] = consensus_df[lfc_cols].std(axis=1)
        
        result = {
            'consensus_genes': consensus_df,
            'de_gene_sets': de_gene_sets,
            'n_per_pipeline': [len(s) for s in de_gene_sets]
        }
        
        return result
    
    def _union(
        self,
        pipeline_results: Dict[str, pd.DataFrame],
        alpha: float,
        lfc_threshold: float
    ) -> Dict:
        """
        Union: genes significant in ANY pipeline.
        
        Liberal approach - includes all genes detected by at least one pipeline.
        """
        logger.info("Using union method (liberal)")
        
        # Get DE genes from each pipeline
        de_gene_sets = []
        all_genes_data = {}
        
        for pipeline, df in pipeline_results.items():
            if df.index.name:
                gene_index = df.index
            else:
                gene_index = df['gene_id']
                df = df.set_index('gene_id')
            
            # Determine DE genes
            if 'padj' in df.columns:
                significant = (df['padj'] < alpha) & (df['log2FoldChange'].abs() > lfc_threshold)
            else:
                significant = (df['pvalue'] < alpha) & (df['log2FoldChange'].abs() > lfc_threshold)
            
            de_genes = set(df.index[significant])
            de_gene_sets.append(de_genes)
            
            # Store data for all DE genes
            for gene in de_genes:
                if gene not in all_genes_data:
                    all_genes_data[gene] = {}
                
                all_genes_data[gene][f'{pipeline}_log2FC'] = df.loc[gene, 'log2FoldChange']
                all_genes_data[gene][f'{pipeline}_pvalue'] = df.loc[gene, 'pvalue']
                all_genes_data[gene][f'{pipeline}_detected'] = True
            
            logger.info(f"{pipeline}: {len(de_genes)} DE genes")
        
        # Find union
        consensus_genes = set.union(*de_gene_sets)
        logger.info(f"Union: {len(consensus_genes)} genes")
        
        # Create consensus DataFrame
        consensus_list = []
        for gene in consensus_genes:
            gene_data = {'gene_id': gene}
            gene_data.update(all_genes_data[gene])
            
            # Count how many pipelines detected this gene
            detected = sum(1 for k in gene_data.keys() if '_detected' in k)
            gene_data['n_pipelines'] = detected
            
            # Calculate mean statistics (only from pipelines that detected it)
            lfc_values = [v for k, v in gene_data.items() if '_log2FC' in k]
            gene_data['mean_log2FC'] = np.mean(lfc_values) if lfc_values else np.nan
            gene_data['std_log2FC'] = np.std(lfc_values) if lfc_values else np.nan
            
            consensus_list.append(gene_data)
        
        consensus_df = pd.DataFrame(consensus_list)
        consensus_df = consensus_df.sort_values('n_pipelines', ascending=False)
        
        result = {
            'consensus_genes': consensus_df,
            'de_gene_sets': de_gene_sets,
            'n_per_pipeline': [len(s) for s in de_gene_sets]
        }
        
        return result
    
    def _calculate_concordance(
        self,
        pipeline_results: Dict[str, pd.DataFrame],
        alpha: float
    ) -> float:
        """Calculate pairwise concordance between pipelines."""
        from itertools import combinations
        
        # Get DE genes from each pipeline
        de_gene_sets = {}
        
        for pipeline, df in pipeline_results.items():
            if df.index.name:
                gene_index = df.index
            else:
                gene_index = df['gene_id']
                df = df.set_index('gene_id')
            
            if 'padj' in df.columns:
                significant = df['padj'] < alpha
            else:
                significant = df['pvalue'] < alpha
            
            de_gene_sets[pipeline] = set(df.index[significant])
        
        # Calculate pairwise Jaccard similarity
        concordances = []
        for pipe1, pipe2 in combinations(de_gene_sets.keys(), 2):
            set1 = de_gene_sets[pipe1]
            set2 = de_gene_sets[pipe2]
            
            if len(set1.union(set2)) > 0:
                jaccard = len(set1.intersection(set2)) / len(set1.union(set2))
                concordances.append(jaccard)
        
        # Return mean concordance
        return np.mean(concordances) if concordances else 0.0
    
    def visualize_concordance(
        self,
        pipeline_results: Dict[str, pd.DataFrame],
        alpha: float = 0.05,
        output_file: Optional[str] = None
    ):
        """
        Visualize concordance between pipelines.
        
        Creates heatmap and Venn diagram showing pipeline agreement.
        """
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        # Get DE genes from each pipeline
        de_gene_sets = {}
        
        for pipeline, df in pipeline_results.items():
            if df.index.name:
                gene_index = df.index
            else:
                gene_index = df['gene_id']
                df = df.set_index('gene_id')
            
            if 'padj' in df.columns:
                significant = df['padj'] < alpha
            else:
                significant = df['pvalue'] < alpha
            
            de_gene_sets[pipeline] = set(df.index[significant])
        
        # Calculate pairwise Jaccard similarity matrix
        pipelines = list(de_gene_sets.keys())
        n_pipelines = len(pipelines)
        
        similarity_matrix = np.zeros((n_pipelines, n_pipelines))
        
        for i, pipe1 in enumerate(pipelines):
            for j, pipe2 in enumerate(pipelines):
                set1 = de_gene_sets[pipe1]
                set2 = de_gene_sets[pipe2]
                
                if i == j:
                    similarity_matrix[i, j] = 1.0
                elif len(set1.union(set2)) > 0:
                    jaccard = len(set1.intersection(set2)) / len(set1.union(set2))
                    similarity_matrix[i, j] = jaccard
        
        # Create heatmap
        fig, ax = plt.subplots(figsize=(10, 8))
        
        sns.heatmap(
            similarity_matrix,
            annot=True,
            fmt='.3f',
            cmap='YlOrRd',
            xticklabels=pipelines,
            yticklabels=pipelines,
            cbar_kws={'label': 'Jaccard Similarity'},
            ax=ax
        )
        
        ax.set_title('Pipeline Concordance (Jaccard Similarity)', fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            logger.info(f"Concordance plot saved to: {output_file}")
        
        plt.close()


if __name__ == '__main__':
    print("RAPTOR Ensemble Analysis")
    print("========================")
    print("\nCombine results from multiple pipelines for robust DE detection.")
    print("\nUsage:")
    print("  from raptor.ensemble_analysis import EnsembleAnalyzer")
    print("  ensemble = EnsembleAnalyzer()")
    print("  results = ensemble.combine_results(pipeline_results, method='weighted_vote')")
