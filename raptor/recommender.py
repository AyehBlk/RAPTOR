"""
Pipeline Recommender

Intelligent recommendation system that matches RNA-seq data characteristics
to optimal analysis pipelines. Uses scoring system to evaluate pipeline
suitability based on data profile.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
"""

import numpy as np
from typing import Dict, List, Optional
import logging

logger = logging.getLogger(__name__)


# Pipeline definitions with characteristics
PIPELINES = {
    1: {
        'name': 'STAR-RSEM-DESeq2',
        'description': 'Gold standard for accuracy',
        'alignment': 'STAR',
        'quantification': 'RSEM',
        'statistics': 'DESeq2',
        'speed': 'slow',
        'memory': 'high',
        'accuracy': 'highest',
        'strengths': [
            'Highest accuracy',
            'Excellent normalization',
            'Handles high variation well',
            'Good for low replication'
        ],
        'best_for': 'Publication-quality results, difficult data'
    },
    2: {
        'name': 'HISAT2-StringTie-Ballgown',
        'description': 'Novel transcript discovery',
        'alignment': 'HISAT2',
        'quantification': 'StringTie',
        'statistics': 'Ballgown',
        'speed': 'medium',
        'memory': 'medium',
        'accuracy': 'medium',
        'strengths': [
            'Novel transcript discovery',
            'Isoform-level analysis',
            'Good for non-model organisms'
        ],
        'best_for': 'Transcriptome assembly, isoform analysis'
    },
    3: {
        'name': 'Salmon-edgeR',
        'description': 'Best balance',
        'alignment': 'Salmon (pseudo)',
        'quantification': 'Salmon',
        'statistics': 'edgeR',
        'speed': 'fast',
        'memory': 'low',
        'accuracy': 'high',
        'strengths': [
            'Excellent balance',
            '3-5× faster than alignment',
            'Low memory usage',
            'Good accuracy'
        ],
        'best_for': 'Most RNA-seq experiments'
    },
    4: {
        'name': 'Kallisto-Sleuth',
        'description': 'Ultra-fast',
        'alignment': 'Kallisto (pseudo)',
        'quantification': 'Kallisto',
        'statistics': 'Sleuth',
        'speed': 'fastest',
        'memory': 'lowest',
        'accuracy': 'good',
        'strengths': [
            'Fastest option',
            'Minimal memory',
            'Good for exploration',
            'Handles large cohorts'
        ],
        'best_for': 'Large datasets, exploratory analysis'
    },
    5: {
        'name': 'STAR-HTSeq-limma-voom',
        'description': 'Flexible modeling',
        'alignment': 'STAR',
        'quantification': 'HTSeq',
        'statistics': 'limma-voom',
        'speed': 'medium',
        'memory': 'high',
        'accuracy': 'high',
        'strengths': [
            'Flexible statistical modeling',
            'Excellent for complex designs',
            'Good batch correction'
        ],
        'best_for': 'Complex experimental designs'
    },
    6: {
        'name': 'STAR-featureCounts-NOISeq',
        'description': 'Non-parametric',
        'alignment': 'STAR',
        'quantification': 'featureCounts',
        'statistics': 'NOISeq',
        'speed': 'slow',
        'memory': 'high',
        'accuracy': 'medium',
        'strengths': [
            'No distribution assumptions',
            'Robust to outliers'
        ],
        'best_for': 'Data not fitting standard distributions'
    },
    7: {
        'name': 'Bowtie2-RSEM-EBSeq',
        'description': 'Bayesian approach',
        'alignment': 'Bowtie2',
        'quantification': 'RSEM',
        'statistics': 'EBSeq',
        'speed': 'very_slow',
        'memory': 'high',
        'accuracy': 'medium',
        'strengths': [
            'Bayesian framework',
            'Handles isoform uncertainty'
        ],
        'best_for': 'Isoform switching analysis'
    },
    8: {
        'name': 'HISAT2-Cufflinks-Cuffdiff',
        'description': 'Legacy comparison',
        'alignment': 'HISAT2',
        'quantification': 'Cufflinks',
        'statistics': 'Cuffdiff',
        'speed': 'slow',
        'memory': 'medium',
        'accuracy': 'low',
        'strengths': [
            'Historical reference',
            'Widely published'
        ],
        'best_for': 'Comparison with legacy studies'
    }
}


class PipelineRecommender:
    """
    Recommend optimal RNA-seq pipeline based on data profile.
    
    Uses a scoring system that considers data characteristics (variation,
    zero-inflation, sample size, depth) and user priorities (accuracy,
    speed, memory) to recommend the best pipeline.
    
    Parameters
    ----------
    profile : dict
        Data profile from RNAseqDataProfiler
    
    Attributes
    ----------
    profile : dict
        Input data profile
    pipelines : dict
        Pipeline definitions
    
    Examples
    --------
    >>> from raptor import RNAseqDataProfiler, PipelineRecommender
    >>> profiler = RNAseqDataProfiler(counts)
    >>> profile = profiler.run_full_profile()
    >>> recommender = PipelineRecommender(profile)
    >>> rec = recommender.get_recommendation(priority='balanced')
    >>> print(rec['primary']['pipeline_name'])
    'Salmon-edgeR'
    """
    
    def __init__(self, profile: Dict):
        """Initialize recommender with data profile."""
        self.profile = profile
        self.pipelines = PIPELINES
        logger.info("Initialized PipelineRecommender")
    
    def get_recommendation(self, priority: str = 'balanced') -> Dict:
        """
        Get pipeline recommendation with reasoning.
        
        Parameters
        ----------
        priority : str
            Optimization priority: 'accuracy', 'speed', 'memory', or 'balanced'
        
        Returns
        -------
        dict
            Recommendation containing:
            - primary: Top recommended pipeline with score and reasoning
            - alternatives: Alternative options
            - all_scores: Scores for all pipelines
        
        Raises
        ------
        ValueError
            If priority is not one of the valid options
        
        Examples
        --------
        >>> rec = recommender.get_recommendation(priority='accuracy')
        >>> print(f"{rec['primary']['pipeline_name']} (Score: {rec['primary']['score']})")
        STAR-RSEM-DESeq2 (Score: 165.3)
        """
        valid_priorities = ['accuracy', 'speed', 'memory', 'balanced']
        if priority not in valid_priorities:
            raise ValueError(f"Priority must be one of {valid_priorities}, got '{priority}'")
        
        logger.info(f"Generating recommendation (priority: {priority})")
        
        # Calculate scores for all pipelines
        all_scores = {}
        for pipeline_id in self.pipelines.keys():
            score, reasoning = self._score_pipeline(pipeline_id, priority)
            all_scores[pipeline_id] = {
                'score': score,
                'reasoning': reasoning,
                'pipeline': self.pipelines[pipeline_id]
            }
        
        # Sort by score
        sorted_pipelines = sorted(all_scores.items(), key=lambda x: x[1]['score'], reverse=True)
        
        # Prepare recommendation
        top_id, top_data = sorted_pipelines[0]
        
        recommendation = {
            'primary': {
                'pipeline_id': top_id,
                'pipeline_name': self.pipelines[top_id]['name'],
                'score': top_data['score'],
                'reasoning': top_data['reasoning'],
                'details': self.pipelines[top_id]
            },
            'alternatives': [],
            'all_scores': {k: v['score'] for k, v in all_scores.items()}
        }
        
        # Add alternatives (top 3)
        for alt_id, alt_data in sorted_pipelines[1:3]:
            recommendation['alternatives'].append({
                'pipeline_id': alt_id,
                'pipeline_name': self.pipelines[alt_id]['name'],
                'score': alt_data['score'],
                'reasoning': alt_data['reasoning'][:2],  # Brief reasoning
                'details': self.pipelines[alt_id]
            })
        
        logger.info(f"Recommended: {recommendation['primary']['pipeline_name']} "
                   f"(Score: {recommendation['primary']['score']:.1f})")
        
        return recommendation
    
    def _score_pipeline(self, pipeline_id: int, priority: str) -> tuple:
        """
        Calculate score for a specific pipeline.
        
        Parameters
        ----------
        pipeline_id : int
            Pipeline ID (1-8)
        priority : str
            User priority
        
        Returns
        -------
        tuple
            (score, reasoning_list)
        
        Notes
        -----
        Scoring components:
        - Data difficulty matching (40%)
        - Sequencing quality matching (30%)
        - Priority weighting (20%)
        - Design complexity (10%)
        
        Score range: 0-200
        """
        pipeline = self.pipelines[pipeline_id]
        score = 0
        reasoning = []
        
        # =====================================================================
        # 1. Data Difficulty Matching (40% weight, max 80 points)
        # =====================================================================
        
        difficulty = self.profile['summary']['difficulty']
        
        if pipeline_id == 1:  # STAR-RSEM-DESeq2
            if difficulty == 'challenging':
                score += 80
                reasoning.append("DESeq2 excels with challenging data")
            elif difficulty == 'moderate':
                score += 60
                reasoning.append("DESeq2 works well with moderate difficulty")
            else:
                score += 40
                reasoning.append("May be overcautious for easy data")
        
        elif pipeline_id == 3:  # Salmon-edgeR
            if difficulty == 'easy':
                score += 80
                reasoning.append("Excellent balance for high-quality data")
            elif difficulty == 'moderate':
                score += 70
                reasoning.append("Good performance with moderate difficulty")
            else:
                score += 50
                reasoning.append("May struggle with very challenging data")
        
        elif pipeline_id == 4:  # Kallisto-Sleuth
            if difficulty == 'easy':
                score += 70
                reasoning.append("Speed advantage with high-quality data")
            elif difficulty == 'moderate':
                score += 50
                reasoning.append("Acceptable for moderate difficulty")
            else:
                score += 20
                reasoning.append("Not recommended for challenging data")
        
        # Library size variation
        cv = self.profile['library_stats']['cv']
        if cv > 0.5:  # High variation
            if pipeline_id in [1, 5]:  # DESeq2, limma-voom
                score += 20
                reasoning.append("Robust normalization for high library variation")
        
        # Zero-inflation
        zero_pct = self.profile['count_distribution']['zero_pct']
        if zero_pct > 60:  # High zero-inflation
            if pipeline_id == 1:  # DESeq2
                score += 15
                reasoning.append("Handles zero-inflation excellently")
        
        # Low replication
        min_reps = self.profile['design']['min_replicates']
        if min_reps < 3:
            if pipeline_id in [1, 7]:  # DESeq2, EBSeq
                score += 15
                reasoning.append("Shrinkage estimators work well with low replication")
        
        # =====================================================================
        # 2. Sequencing Quality Matching (30% weight, max 60 points)
        # =====================================================================
        
        depth = self.profile['sequencing']['depth_category']
        
        if depth == 'low':
            if pipeline_id in [1, 2]:  # Alignment-based
                score += 40
                reasoning.append("Alignment more accurate at low depth")
        elif depth in ['high', 'very_high']:
            if pipeline_id in [3, 4]:  # Pseudo-alignment
                score += 50
                reasoning.append("Pseudo-alignment sufficient at high depth")
            elif pipeline_id == 1:
                score += 30
                reasoning.append("Alignment-based may be overkill at high depth")
        
        # Detection rate
        detection_rate = self.profile['sequencing']['detection_rate']
        if detection_rate > 0.7:
            if pipeline_id in [3, 4]:
                score += 10
                reasoning.append("Good detection enables fast methods")
        
        # =====================================================================
        # 3. Priority Weighting (20% weight, max 40 points)
        # =====================================================================
        
        if priority == 'accuracy':
            accuracy_scores = {'highest': 40, 'high': 30, 'medium': 20, 'good': 25, 'low': 10}
            score += accuracy_scores.get(pipeline['accuracy'], 20)
            if pipeline['accuracy'] == 'highest':
                reasoning.append("Prioritizing accuracy as requested")
        
        elif priority == 'speed':
            speed_scores = {'fastest': 40, 'fast': 35, 'medium': 20, 'slow': 10, 'very_slow': 5}
            score += speed_scores.get(pipeline['speed'], 15)
            if pipeline['speed'] in ['fastest', 'fast']:
                reasoning.append("Fast method matches speed priority")
        
        elif priority == 'memory':
            memory_scores = {'lowest': 40, 'low': 35, 'medium': 20, 'high': 10}
            score += memory_scores.get(pipeline['memory'], 15)
            if pipeline['memory'] in ['lowest', 'low']:
                reasoning.append("Low memory usage matches priority")
        
        else:  # balanced
            # Balanced scoring based on data characteristics
            n_samples = self.profile['design']['n_samples']
            if n_samples > 10 and pipeline_id in [3, 4]:
                score += 30
                reasoning.append("Fast methods appropriate for large sample size")
            elif n_samples < 6 and pipeline_id in [1, 5]:
                score += 35
                reasoning.append("Robust methods important for small samples")
            else:
                score += 20
        
        # =====================================================================
        # 4. Design Complexity (10% weight, max 20 points)
        # =====================================================================
        
        n_conditions = self.profile['design'].get('n_conditions', 2)
        
        if n_conditions > 2:  # Multi-group comparison
            if pipeline_id == 5:  # limma-voom
                score += 20
                reasoning.append("Flexible modeling for complex designs")
            elif pipeline_id in [1, 3]:
                score += 15
            else:
                score += 10
        else:  # Simple two-group
            score += 15  # All methods handle this
        
        # =====================================================================
        # Penalties
        # =====================================================================
        
        # Penalize outdated methods
        if pipeline_id == 8:  # Cuffdiff
            score -= 20
            reasoning.append("Note: This is a legacy method")
        
        # Penalize very slow methods for large datasets
        if n_samples > 20 and pipeline['speed'] in ['slow', 'very_slow']:
            score -= 15
            reasoning.append("Consider faster alternatives for large dataset")
        
        # Ensure score is in range [0, 200]
        score = max(0, min(200, score))
        
        return score, reasoning
    
    def compare_pipelines(self, pipeline_ids: List[int], priority: str = 'balanced') -> Dict:
        """
        Compare specific pipelines side-by-side.
        
        Parameters
        ----------
        pipeline_ids : list of int
            Pipeline IDs to compare
        priority : str
            Optimization priority
        
        Returns
        -------
        dict
            Comparison with scores and reasoning for each pipeline
        
        Examples
        --------
        >>> comparison = recommender.compare_pipelines([1, 3, 4])
        >>> for p_id, data in comparison.items():
        ...     print(f"Pipeline {p_id}: {data['score']}")
        Pipeline 1: 145.5
        Pipeline 3: 162.3
        Pipeline 4: 125.7
        """
        comparison = {}
        
        for pipeline_id in pipeline_ids:
            if pipeline_id not in self.pipelines:
                logger.warning(f"Invalid pipeline ID: {pipeline_id}")
                continue
            
            score, reasoning = self._score_pipeline(pipeline_id, priority)
            comparison[pipeline_id] = {
                'pipeline': self.pipelines[pipeline_id],
                'score': score,
                'reasoning': reasoning
            }
        
        return comparison


# =============================================================================
# Convenience Functions
# =============================================================================

def quick_recommend(profile: Dict, priority: str = 'balanced') -> str:
    """
    Get quick recommendation as pipeline name.
    
    Parameters
    ----------
    profile : dict
        Data profile
    priority : str
        Optimization priority
    
    Returns
    -------
    str
        Recommended pipeline name
    
    Examples
    --------
    >>> pipeline_name = quick_recommend(profile, priority='speed')
    >>> print(pipeline_name)
    'Kallisto-Sleuth'
    """
    recommender = PipelineRecommender(profile)
    rec = recommender.get_recommendation(priority)
    return rec['primary']['pipeline_name']


if __name__ == '__main__':
    print("RAPTOR Pipeline Recommender")
    print("===========================")
    print("\nAvailable pipelines:")
    for pid, pipeline in PIPELINES.items():
        print(f"  {pid}. {pipeline['name']} - {pipeline['description']}")
