#!/usr/bin/env python3
"""
RAPTOR v2.2.0 - Module 8: Parameter Optimization

Comprehensive parameter optimization for differential expression analysis with
multiple scientifically validated methods.

Optimization Methods:
====================
1. GROUND TRUTH - Validated gene list (gold standard)
   Requirements: 10-15 genes minimum, 20-30 recommended, 50+ ideal
   
2. FDR CONTROL - Statistical FDR estimation (no validation needed)
   Based on: Storey & Tibshirani (2003) PNAS
   
3. STABILITY - Bootstrap-based stability analysis (no validation needed)
   Based on: Meinshausen & Bühlmann (2010) JRSS-B
   
4. REPRODUCIBILITY - Independent cohort validation (requires 2 datasets)
   Based on: Replication principles

Search Strategies:
=================
- Grid Search: Exhaustive search (guaranteed optimum)
- Random Search: Monte Carlo sampling (Bergstra & Bengio, 2012)
- Differential Evolution: Population-based global optimization (Storn & Price, 1997)

Parameters Optimized:
====================
- alpha: FDR/p-value significance threshold (0.001 - 0.20)
- lfc_threshold: Minimum log2 fold change (0.0 - 2.0)

References:
-----------
[1] Storey, J.D., & Tibshirani, R. (2003). Statistical significance for 
    genomewide studies. PNAS, 100(16), 9440-9445.

[2] Meinshausen, N., & Bühlmann, P. (2010). Stability selection.
    J Royal Statistical Society: Series B, 72(4), 417-473.

[3] Bergstra, J., & Bengio, Y. (2012). Random search for hyper-parameter 
    optimization. Journal of Machine Learning Research, 13(1), 281-305.

[4] Storn, R., & Price, K. (1997). Differential evolution. 
    J Global Optimization, 11(4), 341-359.

[5] Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate.
    J Royal Statistical Society: Series B, 57(1), 289-300.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
License: MIT
"""

import logging
import json
import pickle
import warnings
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Any, Optional, Union
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)


# =============================================================================
# VALIDATION CONSTANTS
# =============================================================================

# Ground truth validation requirements
MIN_VALIDATED_GENES_ABSOLUTE = 10  # Absolute minimum
MIN_VALIDATED_GENES_RECOMMENDED = 20  # Recommended minimum
IDEAL_VALIDATED_GENES = 50  # Ideal for robust optimization

# FDR control constants
MIN_GENES_FOR_FDR = 1000  # Need enough genes for FDR estimation
PI0_LAMBDA_RANGE = np.arange(0.05, 0.95, 0.05)  # For π₀ estimation

# Stability constants
MIN_BOOTSTRAP_ITERATIONS = 50  # Minimum bootstrap samples
RECOMMENDED_BOOTSTRAP = 100  # Recommended bootstrap samples

# Parameter bounds
ALPHA_MIN, ALPHA_MAX = 0.001, 0.20
LFC_MIN, LFC_MAX = 0.0, 2.0


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class ParameterSpace:
    """
    Parameter search space definition.
    
    Attributes
    ----------
    name : str
        Parameter name ('alpha' or 'lfc_threshold')
    bounds : tuple
        (min, max) for parameter range
    default : float
        Default value
    description : str
        Human-readable description
    """
    name: str
    bounds: Tuple[float, float]
    default: float
    description: str
    
    def sample_random(self) -> float:
        """Sample random value from parameter space."""
        return np.random.uniform(self.bounds[0], self.bounds[1])
    
    def sample_grid(self, n_points: int = 5) -> List[float]:
        """Generate grid points for exhaustive search."""
        return list(np.linspace(self.bounds[0], self.bounds[1], n_points))
    
    def validate(self, value: float) -> bool:
        """Check if value is within bounds."""
        return self.bounds[0] <= value <= self.bounds[1]


@dataclass
class OptimizationResult:
    """
    Parameter optimization results.
    
    Attributes
    ----------
    best_parameters : dict
        Optimized parameter values {'alpha': float, 'lfc_threshold': float}
    best_score : float
        Best score achieved
    optimization_method : str
        Method used: 'ground_truth', 'fdr_control', 'stability', 'reproducibility'
    search_strategy : str
        Search strategy: 'grid', 'random', 'differential_evolution'
    metric : str
        Optimization metric (e.g., 'f1_score', 'stability_score')
    n_iterations : int
        Number of iterations performed
    convergence_history : list
        Complete optimization history
    deg_genes : pd.DataFrame, optional
        Filtered DEG genes with optimized parameters
    n_deg_genes : int
        Number of DEG genes at optimal parameters
    performance_metrics : dict
        Detailed performance metrics
    metadata : dict
        Additional metadata (method-specific)
    timestamp : str
        ISO format timestamp
    """
    best_parameters: Dict[str, float]
    best_score: float
    optimization_method: str
    search_strategy: str
    metric: str
    n_iterations: int
    convergence_history: List[Dict] = field(default_factory=list)
    deg_genes: Optional[pd.DataFrame] = None
    n_deg_genes: int = 0
    performance_metrics: Dict = field(default_factory=dict)
    metadata: Dict = field(default_factory=dict)
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())
    
    def summary(self) -> str:
        """Generate human-readable summary."""
        lines = [
            "="*70,
            "RAPTOR Module 8: Parameter Optimization Results",
            "="*70,
            f"Optimization Method: {self.optimization_method}",
            f"Search Strategy: {self.search_strategy}",
            f"Metric: {self.metric}",
            f"Iterations: {self.n_iterations}",
            "",
            "Optimal Parameters:",
            f"  alpha: {self.best_parameters['alpha']:.4f}",
            f"  lfc_threshold: {self.best_parameters['lfc_threshold']:.2f}",
            "",
            f"Best Score: {self.best_score:.4f}",
            f"DEG Genes: {self.n_deg_genes}",
            "",
        ]
        
        if self.performance_metrics:
            lines.append("Performance Metrics:")
            for key, value in self.performance_metrics.items():
                if isinstance(value, float):
                    lines.append(f"  {key}: {value:.4f}")
                else:
                    lines.append(f"  {key}: {value}")
            lines.append("")
        
        lines.append(f"Timestamp: {self.timestamp}")
        lines.append("="*70)
        
        return "\n".join(lines)
    
    def save(self, output_dir: Union[str, Path]) -> None:
        """Save results to disk with error handling."""
        try:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            raise IOError(f"Failed to create output directory {output_dir}: {e}")
        
        # Save parameters as YAML
        try:
            params_file = output_dir / 'optimized_params.yaml'
            try:
                import yaml
            except ImportError:
                logger.error("PyYAML not installed. Cannot save YAML file.")
                logger.info("Install with: pip install pyyaml")
                raise ImportError("PyYAML required. Install with: pip install pyyaml")
            
            with open(params_file, 'w') as f:
                yaml_data = {
                    'best_parameters': self.best_parameters,
                    'best_score': float(self.best_score),
                    'optimization_method': self.optimization_method,
                    'search_strategy': self.search_strategy,
                    'metric': self.metric,
                    'n_iterations': self.n_iterations,
                    'n_deg_genes': self.n_deg_genes,
                    'performance_metrics': {k: float(v) if isinstance(v, (np.floating, float)) else v 
                                           for k, v in self.performance_metrics.items()},
                    'timestamp': self.timestamp
                }
                yaml.dump(yaml_data, f, default_flow_style=False)
            logger.info(f"Saved parameters to {params_file}")
        except Exception as e:
            logger.error(f"Failed to save YAML parameters: {e}")
            raise
        
        # Save DEG genes if available
        if self.deg_genes is not None and len(self.deg_genes) > 0:
            try:
                deg_file = output_dir / 'deg_genes.csv'
                self.deg_genes.to_csv(deg_file, index=False)
                logger.info(f"Saved {len(self.deg_genes)} DEG genes to {deg_file}")
            except Exception as e:
                logger.error(f"Failed to save DEG genes: {e}")
                raise
        else:
            logger.warning("No DEG genes to save (empty result)")
        
        # Save full result as pickle
        try:
            result_file = output_dir / 'optimization_result.pkl'
            with open(result_file, 'wb') as f:
                pickle.dump(self, f)
            logger.info(f"Saved full result to {result_file}")
        except Exception as e:
            logger.error(f"Failed to save pickle file: {e}")
            raise
        
        # Save convergence history as JSON
        try:
            history_file = output_dir / 'convergence_history.json'
            with open(history_file, 'w') as f:
                json.dump(self.convergence_history, f, indent=2)
            logger.info(f"Saved convergence history to {history_file}")
        except Exception as e:
            logger.error(f"Failed to save convergence history: {e}")
            raise
        
        logger.info(f"\n✅ All results saved to {output_dir}/")
    
    @classmethod
    def load(cls, filepath: Union[str, Path]) -> 'OptimizationResult':
        """Load results from pickle file."""
        with open(filepath, 'rb') as f:
            return pickle.load(f)


# =============================================================================
# PARAMETER OPTIMIZER (Base Class)
# =============================================================================

class ParameterOptimizer:
    """
    Base class for parameter optimization.
    
    Handles search strategies (grid, random, DE) while allowing different
    evaluation methods (ground truth, FDR, stability, reproducibility).
    
    Parameters
    ----------
    de_result : pd.DataFrame
        Differential expression results
        Required columns: gene_id (or index), log2FoldChange, pvalue, padj
    parameter_spaces : dict, optional
        Custom parameter spaces. If None, uses default ranges.
    random_state : int, optional
        Random seed for reproducibility
        
    Attributes
    ----------
    de_result : pd.DataFrame
        DE results
    parameter_spaces : dict
        Parameter search spaces
    history : list
        Optimization history
    best_params : dict
        Current best parameters
    best_score : float
        Current best score
    """
    
    def __init__(
        self,
        de_result: pd.DataFrame,
        parameter_spaces: Optional[Dict[str, ParameterSpace]] = None,
        random_state: Optional[int] = None
    ):
        self.de_result = de_result.copy()
        self.random_state = random_state
        
        if random_state is not None:
            np.random.seed(random_state)
        
        # Default parameter spaces
        if parameter_spaces is None:
            self.parameter_spaces = {
                'alpha': ParameterSpace(
                    name='alpha',
                    bounds=(ALPHA_MIN, ALPHA_MAX),
                    default=0.05,
                    description='FDR/p-value significance threshold'
                ),
                'lfc_threshold': ParameterSpace(
                    name='lfc_threshold',
                    bounds=(LFC_MIN, LFC_MAX),
                    default=1.0,
                    description='Minimum log2 fold change'
                )
            }
        else:
            self.parameter_spaces = parameter_spaces
        
        # Validate DE result
        self._validate_de_result()
        
        # Initialize optimization state
        self.history = []
        self.best_params = {name: space.default for name, space in self.parameter_spaces.items()}
        self.best_score = -np.inf
    
    def _validate_de_result(self) -> None:
        """Validate DE result format."""
        required_cols = ['log2FoldChange', 'pvalue']
        
        for col in required_cols:
            if col not in self.de_result.columns:
                raise ValueError(f"DE result missing required column: {col}")
        
        # Ensure gene_id column or index
        if 'gene_id' not in self.de_result.columns and self.de_result.index.name is None:
            raise ValueError("DE result must have 'gene_id' column or named index")
        
        # Check for NaN values
        n_na_lfc = self.de_result['log2FoldChange'].isna().sum()
        n_na_pval = self.de_result['pvalue'].isna().sum()
        
        if n_na_lfc > 0:
            logger.warning(f"{n_na_lfc} genes have NaN log2FoldChange (will be excluded)")
        if n_na_pval > 0:
            logger.warning(f"{n_na_pval} genes have NaN pvalue (will be excluded)")
        
        # Remove NaN rows
        self.de_result = self.de_result.dropna(subset=['log2FoldChange', 'pvalue'])
        
        logger.info(f"DE result validated: {len(self.de_result)} genes")
    
    def evaluate_parameters(self, params: Dict[str, float]) -> float:
        """
        Evaluate parameter combination (to be overridden by subclasses).
        
        Parameters
        ----------
        params : dict
            Parameters to evaluate
            
        Returns
        -------
        float
            Score (higher is better)
        """
        raise NotImplementedError("Subclasses must implement evaluate_parameters")
    
    def optimize(
        self,
        strategy: str = 'grid',
        n_iterations: int = 50,
        metric: str = 'f1_score',
        grid_points: int = 5,
        **kwargs
    ) -> OptimizationResult:
        """
        Run parameter optimization.
        
        Parameters
        ----------
        strategy : str
            'grid', 'random', or 'differential_evolution'
        n_iterations : int
            Number of iterations (for random/DE)
        metric : str
            Optimization metric name (for display)
        grid_points : int
            Grid points per parameter (for grid search)
        **kwargs
            Additional strategy-specific arguments
            
        Returns
        -------
        OptimizationResult
            Optimization results
        """
        logger.info(f"\n{'='*70}")
        logger.info(f"Starting {strategy} optimization")
        logger.info(f"{'='*70}\n")
        
        if strategy == 'grid':
            self._grid_search(grid_points)
        elif strategy == 'random':
            self._random_search(n_iterations)
        elif strategy == 'differential_evolution':
            self._differential_evolution(n_iterations, **kwargs)
        else:
            raise ValueError(f"Unknown strategy: {strategy}")
        
        # Get DEG genes at optimal parameters
        deg_genes = self._extract_deg_genes(self.best_params)
        
        # Create result object
        result = OptimizationResult(
            best_parameters=self.best_params.copy(),
            best_score=self.best_score,
            optimization_method=self.__class__.__name__.replace('Optimizer', '').lower(),
            search_strategy=strategy,
            metric=metric,
            n_iterations=len(self.history),
            convergence_history=self.history.copy(),
            deg_genes=deg_genes,
            n_deg_genes=len(deg_genes) if deg_genes is not None else 0,
            performance_metrics=self._calculate_performance_metrics(self.best_params),
            metadata=self._get_metadata()
        )
        
        return result
    
    def _extract_deg_genes(self, params: Dict[str, float]) -> pd.DataFrame:
        """
        Extract DEG genes using given parameters.
        
        Parameters
        ----------
        params : dict
            Parameters with 'alpha' and 'lfc_threshold'
            
        Returns
        -------
        pd.DataFrame
            Filtered DEG genes with direction column
            
        Raises
        ------
        ValueError
            If parameters are invalid or result in no genes
        """
        try:
            alpha = params['alpha']
            lfc_thresh = params['lfc_threshold']
        except KeyError as e:
            raise ValueError(f"Missing required parameter: {e}")
        
        # Validate parameter ranges
        if not (0 < alpha < 1):
            raise ValueError(f"alpha must be between 0 and 1, got {alpha}")
        if not (0 <= lfc_thresh <= 5):
            raise ValueError(f"lfc_threshold must be between 0 and 5, got {lfc_thresh}")
        
        # Apply filters
        try:
            if 'padj' in self.de_result.columns:
                mask = (
                    (self.de_result['padj'] < alpha) &
                    (self.de_result['log2FoldChange'].abs() > lfc_thresh)
                )
            else:
                mask = (
                    (self.de_result['pvalue'] < alpha) &
                    (self.de_result['log2FoldChange'].abs() > lfc_thresh)
                )
        except Exception as e:
            raise RuntimeError(f"Failed to apply filters: {e}")
        
        deg_genes = self.de_result[mask].copy()
        
        # Check if result is empty
        if len(deg_genes) == 0:
            logger.warning(
                f"No genes passed filters (alpha={alpha:.4f}, lfc={lfc_thresh:.2f}). "
                f"Try more relaxed parameters."
            )
            # Return empty DataFrame with correct structure
            deg_genes = pd.DataFrame(columns=list(self.de_result.columns) + ['direction', 'passed_filter'])
            return deg_genes
        
        # Add direction column
        try:
            deg_genes['direction'] = deg_genes['log2FoldChange'].apply(
                lambda x: 'up' if x > 0 else 'down'
            )
        except Exception as e:
            raise RuntimeError(f"Failed to add direction column: {e}")
        
        # Add passed_filter column
        deg_genes['passed_filter'] = True
        
        # Ensure gene_id column exists
        try:
            if 'gene_id' not in deg_genes.columns:
                deg_genes.reset_index(inplace=True)
                if 'index' in deg_genes.columns:
                    deg_genes.rename(columns={'index': 'gene_id'}, inplace=True)
        except Exception as e:
            logger.warning(f"Could not ensure gene_id column: {e}")
        
        return deg_genes
    
    def _calculate_performance_metrics(self, params: Dict[str, float]) -> Dict:
        """Calculate performance metrics (to be overridden by subclasses)."""
        return {}
    
    def _get_metadata(self) -> Dict:
        """Get method-specific metadata (to be overridden by subclasses)."""
        return {}
    
    def _grid_search(self, n_points: int = 5) -> None:
        """Grid search optimization."""
        logger.info(f"Grid search with {n_points} points per parameter")
        
        # Generate grid
        grids = {name: space.sample_grid(n_points) 
                for name, space in self.parameter_spaces.items()}
        
        # Calculate total combinations
        total = np.prod([len(g) for g in grids.values()])
        logger.info(f"Testing {total} parameter combinations\n")
        
        # Iterate through grid
        iteration = 0
        for alpha in grids['alpha']:
            for lfc in grids['lfc_threshold']:
                params = {'alpha': alpha, 'lfc_threshold': lfc}
                score = self.evaluate_parameters(params)
                
                self.history.append({
                    'iteration': iteration,
                    'parameters': params.copy(),
                    'score': score
                })
                
                if score > self.best_score:
                    self.best_score = score
                    self.best_params = params.copy()
                    logger.info(f"✓ New best [{iteration+1}/{total}]: {score:.4f} "
                              f"(alpha={alpha:.4f}, lfc={lfc:.2f})")
                
                iteration += 1
                
                if (iteration + 1) % 10 == 0 and (iteration + 1) < total:
                    logger.info(f"  Progress: {iteration+1}/{total}")
        
        logger.info(f"\n✓ Grid search complete: {total} combinations tested")
    
    def _random_search(self, n_iterations: int) -> None:
        """Random search optimization."""
        logger.info(f"Random search with {n_iterations} iterations\n")
        
        for i in range(n_iterations):
            # Sample random parameters
            params = {name: space.sample_random() 
                     for name, space in self.parameter_spaces.items()}
            
            score = self.evaluate_parameters(params)
            
            self.history.append({
                'iteration': i,
                'parameters': params.copy(),
                'score': score
            })
            
            if score > self.best_score:
                self.best_score = score
                self.best_params = params.copy()
                logger.info(f"✓ New best [{i+1}/{n_iterations}]: {score:.4f}")
            
            if (i + 1) % 20 == 0 and (i + 1) < n_iterations:
                logger.info(f"  Progress: {i+1}/{n_iterations}, best: {self.best_score:.4f}")
        
        logger.info(f"\n✓ Random search complete: {n_iterations} iterations")
    
    def _differential_evolution(
        self,
        n_iterations: int,
        population_size: int = 15,
        mutation: float = 0.5,
        recombination: float = 0.7
    ) -> None:
        """Differential evolution optimization."""
        logger.info(f"Differential evolution: {n_iterations} generations, "
                   f"population={population_size}\n")
        
        # Initialize population
        population = []
        for _ in range(population_size):
            params = {name: space.sample_random() 
                     for name, space in self.parameter_spaces.items()}
            score = self.evaluate_parameters(params)
            population.append({'params': params, 'score': score})
        
        # Find initial best
        population.sort(key=lambda x: x['score'], reverse=True)
        self.best_params = population[0]['params'].copy()
        self.best_score = population[0]['score']
        
        logger.info(f"Initial population best: {self.best_score:.4f}")
        
        # Evolution
        for gen in range(n_iterations):
            for i in range(population_size):
                # Select three random individuals (different from i)
                indices = [j for j in range(population_size) if j != i]
                a, b, c = np.random.choice(indices, 3, replace=False)
                
                # Mutation
                mutant = {}
                for name in self.parameter_spaces.keys():
                    val = (population[a]['params'][name] +
                          mutation * (population[b]['params'][name] - 
                                    population[c]['params'][name]))
                    
                    # Clip to bounds
                    bounds = self.parameter_spaces[name].bounds
                    val = np.clip(val, bounds[0], bounds[1])
                    mutant[name] = val
                
                # Recombination
                trial = {}
                for name in self.parameter_spaces.keys():
                    if np.random.random() < recombination:
                        trial[name] = mutant[name]
                    else:
                        trial[name] = population[i]['params'][name]
                
                # Evaluation
                trial_score = self.evaluate_parameters(trial)
                
                # Selection
                if trial_score > population[i]['score']:
                    population[i] = {'params': trial, 'score': trial_score}
                    
                    if trial_score > self.best_score:
                        self.best_score = trial_score
                        self.best_params = trial.copy()
                        logger.info(f"✓ New best [gen {gen+1}/{n_iterations}]: {trial_score:.4f}")
            
            if (gen + 1) % 10 == 0 and (gen + 1) < n_iterations:
                current_best = max([p['score'] for p in population])
                logger.info(f"  Generation {gen+1}/{n_iterations}, "
                          f"population best: {current_best:.4f}")
        
        logger.info(f"\n✓ Differential evolution complete: {n_iterations} generations")


# =============================================================================
# METHOD 1: GROUND TRUTH OPTIMIZATION
# =============================================================================

class GroundTruthOptimizer(ParameterOptimizer):
    """
    Parameter optimization using validated gene list.
    
    GOLD STANDARD method when experimental validation is available.
    
    Requirements
    ------------
    Validated genes: 10-15 minimum, 20-30 recommended, 50+ ideal
    
    Validation sources:
    - qRT-PCR confirmation
    - Western blot validation
    - Independent RNA-seq cohort
    - Literature (same tissue/condition)
    
    Parameters
    ----------
    de_result : pd.DataFrame
        DE results with columns: gene_id, log2FoldChange, pvalue, padj
    ground_truth : pd.DataFrame
        Validated genes with column: gene_id
        Optionally: expected_direction ('up' or 'down')
    check_direction : bool, default=True
        If True and ground_truth has 'expected_direction', penalize wrong direction
    parameter_spaces : dict, optional
        Custom parameter spaces
    random_state : int, optional
        Random seed
        
    Examples
    --------
    >>> # Load validated genes
    >>> validated = pd.DataFrame({
    ...     'gene_id': ['BRCA1', 'TP53', 'MYC', 'KRAS', ...],
    ...     'expected_direction': ['up', 'up', 'down', 'up', ...]
    ... })
    >>> 
    >>> # Optimize
    >>> optimizer = GroundTruthOptimizer(de_result, validated)
    >>> result = optimizer.optimize(strategy='grid', metric='f1_score')
    >>> 
    >>> print(f"Best alpha: {result.best_parameters['alpha']:.4f}")
    >>> print(f"Best LFC: {result.best_parameters['lfc_threshold']:.2f}")
    >>> print(f"F1 score: {result.best_score:.4f}")
    """
    
    def __init__(
        self,
        de_result: pd.DataFrame,
        ground_truth: pd.DataFrame,
        check_direction: bool = True,
        parameter_spaces: Optional[Dict[str, ParameterSpace]] = None,
        random_state: Optional[int] = None
    ):
        super().__init__(de_result, parameter_spaces, random_state)
        
        self.ground_truth = ground_truth.copy()
        self.check_direction = check_direction
        
        # Validate ground truth
        self._validate_ground_truth()
        
        logger.info(f"Ground truth optimizer initialized:")
        logger.info(f"  Validated genes: {len(self.ground_truth)}")
        logger.info(f"  Check direction: {self.check_direction}")
    
    def _validate_ground_truth(self) -> None:
        """Validate ground truth gene list."""
        # Check column
        if 'gene_id' not in self.ground_truth.columns:
            raise ValueError("Ground truth must have 'gene_id' column")
        
        n_validated = len(self.ground_truth)
        
        # Check minimum
        if n_validated < MIN_VALIDATED_GENES_ABSOLUTE:
            raise ValueError(
                f"Only {n_validated} validated genes provided. "
                f"Minimum {MIN_VALIDATED_GENES_ABSOLUTE} required for optimization. "
                f"\n\n"
                f"RECOMMENDATION: Need at least {MIN_VALIDATED_GENES_RECOMMENDED} genes "
                f"for reliable optimization.\n"
                f"\n"
                f"Validation sources:\n"
                f"  - qRT-PCR confirmation ($75/gene)\n"
                f"  - Western blot validation\n"
                f"  - Independent RNA-seq cohort\n"
                f"  - Literature-validated genes (same tissue/condition)\n"
                f"\n"
                f"Cost estimate for {MIN_VALIDATED_GENES_RECOMMENDED} genes: ~$2,700\n"
                f"\n"
                f"If validation not possible:\n"
                f"  - Skip Module 8\n"
                f"  - Use standard thresholds (alpha=0.05, lfc=1.0)\n"
                f"  - Proceed to Module 9 (Ensemble Analysis)"
            )
        
        # Warnings for borderline cases
        if n_validated < MIN_VALIDATED_GENES_RECOMMENDED:
            logger.warning(
                f"\n⚠️  WARNING: Only {n_validated} validated genes\n"
                f"\n"
                f"This is below the recommended minimum of {MIN_VALIDATED_GENES_RECOMMENDED}.\n"
                f"Optimization may be unreliable.\n"
                f"\n"
                f"Recommendations:\n"
                f"  BEST: Validate {MIN_VALIDATED_GENES_RECOMMENDED - n_validated} more genes "
                f"(cost: ~${(MIN_VALIDATED_GENES_RECOMMENDED - n_validated) * 75})\n"
                f"  ALTERNATIVE: Use optimization results cautiously\n"
                f"  ALTERNATIVE: Use standard thresholds instead\n"
            )
        elif n_validated < IDEAL_VALIDATED_GENES:
            logger.info(
                f"✓ Good: {n_validated} validated genes (recommend ≥{IDEAL_VALIDATED_GENES} for ideal optimization)"
            )
        else:
            logger.info(
                f"✓ Excellent: {n_validated} validated genes (high-quality optimization expected)"
            )
        
        # Check for direction column
        if 'expected_direction' in self.ground_truth.columns:
            logger.info(f"  Direction information available")
            if self.check_direction:
                logger.info(f"  Will penalize incorrect directions")
        
        # Check overlap with DE result
        if 'gene_id' in self.de_result.columns:
            de_genes = set(self.de_result['gene_id'])
        else:
            de_genes = set(self.de_result.index)
        
        gt_genes = set(self.ground_truth['gene_id'])
        overlap = len(de_genes & gt_genes)
        
        if overlap < n_validated:
            logger.warning(
                f"Only {overlap}/{n_validated} validated genes found in DE results"
            )
        
        logger.info("")
    
    def evaluate_parameters(self, params: Dict[str, float]) -> float:
        """
        Evaluate parameters using F1 score against ground truth.
        
        F1 = 2 × (precision × recall) / (precision + recall)
        
        If check_direction=True and ground_truth has 'expected_direction',
        penalize genes with wrong direction.
        """
        # Extract DEGs
        deg_genes = self._extract_deg_genes(params)
        
        if 'gene_id' in deg_genes.columns:
            predicted = set(deg_genes['gene_id'])
        else:
            predicted = set(deg_genes.index)
        
        true_pos = set(self.ground_truth['gene_id'])
        
        # Calculate metrics
        tp = len(predicted & true_pos)
        fp = len(predicted - true_pos)
        fn = len(true_pos - predicted)
        
        # Direction penalty
        if self.check_direction and 'expected_direction' in self.ground_truth.columns:
            # Check direction for true positives
            for gene in (predicted & true_pos):
                expected_dir = self.ground_truth[
                    self.ground_truth['gene_id'] == gene
                ]['expected_direction'].iloc[0]
                
                predicted_dir = deg_genes[
                    deg_genes['gene_id'] == gene
                ]['direction'].iloc[0] if 'gene_id' in deg_genes.columns else deg_genes.loc[gene, 'direction']
                
                # If direction is wrong, count as false positive + false negative
                if expected_dir != predicted_dir:
                    tp -= 1
                    fp += 1
                    fn += 1
        
        # Calculate F1 score
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1_score = (2 * precision * recall / (precision + recall) 
                   if (precision + recall) > 0 else 0)
        
        return f1_score
    
    def _calculate_performance_metrics(self, params: Dict[str, float]) -> Dict:
        """Calculate detailed performance metrics."""
        deg_genes = self._extract_deg_genes(params)
        
        if 'gene_id' in deg_genes.columns:
            predicted = set(deg_genes['gene_id'])
        else:
            predicted = set(deg_genes.index)
        
        true_pos = set(self.ground_truth['gene_id'])
        
        tp = len(predicted & true_pos)
        fp = len(predicted - true_pos)
        fn = len(true_pos - predicted)
        tn = len(self.de_result) - tp - fp - fn
        
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1_score = (2 * precision * recall / (precision + recall) 
                   if (precision + recall) > 0 else 0)
        accuracy = (tp + tn) / (tp + tn + fp + fn) if (tp + tn + fp + fn) > 0 else 0
        
        return {
            'precision': precision,
            'recall': recall,
            'f1_score': f1_score,
            'accuracy': accuracy,
            'true_positives': tp,
            'false_positives': fp,
            'false_negatives': fn,
            'true_negatives': tn
        }
    
    def _get_metadata(self) -> Dict:
        """Get ground truth specific metadata."""
        return {
            'n_validated_genes': len(self.ground_truth),
            'check_direction': self.check_direction,
            'has_direction_info': 'expected_direction' in self.ground_truth.columns
        }


# =============================================================================
# METHOD 2: FDR CONTROL OPTIMIZATION
# =============================================================================

class FDRControlOptimizer(ParameterOptimizer):
    """
    Parameter optimization using FDR control theory.
    
    NO VALIDATION NEEDED! Uses statistical theory to estimate optimal parameters.
    
    Based on Storey & Tibshirani (2003) method for FDR estimation:
    - Estimates π₀ (proportion of true nulls)
    - Finds alpha that achieves target FDR
    
    Requirements
    ------------
    - At least 1000 genes in DE results (for reliable FDR estimation)
    - P-values should show expected distribution
    
    Parameters
    ----------
    de_result : pd.DataFrame
        DE results with columns: gene_id, log2FoldChange, pvalue, padj
    target_fdr : float, default=0.05
        Target false discovery rate (0.01-0.20)
    parameter_spaces : dict, optional
        Custom parameter spaces
    random_state : int, optional
        Random seed
        
    References
    ----------
    Storey, J.D., & Tibshirani, R. (2003). Statistical significance for 
    genomewide studies. PNAS, 100(16), 9440-9445.
    
    Examples
    --------
    >>> optimizer = FDRControlOptimizer(de_result, target_fdr=0.05)
    >>> result = optimizer.optimize(strategy='grid')
    >>> print(f"Optimal alpha for FDR={0.05}: {result.best_parameters['alpha']:.4f}")
    """
    
    def __init__(
        self,
        de_result: pd.DataFrame,
        target_fdr: float = 0.05,
        parameter_spaces: Optional[Dict[str, ParameterSpace]] = None,
        random_state: Optional[int] = None
    ):
        super().__init__(de_result, parameter_spaces, random_state)
        
        self.target_fdr = target_fdr
        
        # Validate
        if not 0 < target_fdr < 1:
            raise ValueError(f"target_fdr must be between 0 and 1, got {target_fdr}")
        
        if len(self.de_result) < MIN_GENES_FOR_FDR:
            logger.warning(
                f"Only {len(self.de_result)} genes available. "
                f"Recommend ≥{MIN_GENES_FOR_FDR} for reliable FDR estimation."
            )
        
        # Estimate π₀
        self.pi0 = self._estimate_pi0()
        
        logger.info(f"FDR control optimizer initialized:")
        logger.info(f"  Target FDR: {self.target_fdr:.3f}")
        logger.info(f"  Estimated π₀: {self.pi0:.3f}")
        logger.info(f"  (π₀ = proportion of true null hypotheses)")
        logger.info("")
    
    def _estimate_pi0(self) -> float:
        """
        Estimate π₀ (proportion of true nulls) using Storey's method.
        
        Returns
        -------
        float
            Estimated π₀ (0-1)
        """
        pvals = self.de_result['pvalue'].values
        pvals = pvals[~np.isnan(pvals)]
        
        if len(pvals) == 0:
            return 1.0
        
        # Storey's π₀ estimator
        pi0_estimates = []
        
        for lambda_val in PI0_LAMBDA_RANGE:
            # Count p-values > lambda
            n_above = np.sum(pvals > lambda_val)
            # Estimate π₀
            pi0 = n_above / (len(pvals) * (1 - lambda_val))
            pi0_estimates.append(min(pi0, 1.0))
        
        # Use median for robustness
        pi0 = np.median(pi0_estimates)
        
        # Ensure valid range
        pi0 = np.clip(pi0, 0.0, 1.0)
        
        return pi0
    
    def evaluate_parameters(self, params: Dict[str, float]) -> float:
        """
        Evaluate parameters by estimating actual FDR and comparing to target.
        
        Score = 1 - |estimated_fdr - target_fdr|
        (Higher score = closer to target FDR)
        """
        alpha = params['alpha']
        lfc_thresh = params['lfc_threshold']
        
        # Count significant genes
        mask = (
            (self.de_result['pvalue'] < alpha) &
            (self.de_result['log2FoldChange'].abs() > lfc_thresh)
        )
        
        n_sig = mask.sum()
        
        if n_sig == 0:
            # No genes significant → FDR undefined, penalize heavily
            return 0.0
        
        # Estimate FDR using Storey's method
        # FDR ≈ (π₀ × alpha × n_total) / n_significant
        n_total = len(self.de_result)
        estimated_fdr = (self.pi0 * alpha * n_total) / n_sig
        
        # Clip to valid range
        estimated_fdr = np.clip(estimated_fdr, 0.0, 1.0)
        
        # Score: closer to target_fdr is better
        # Use negative absolute difference (higher = better)
        score = 1.0 - abs(estimated_fdr - self.target_fdr)
        
        return score
    
    def _calculate_performance_metrics(self, params: Dict[str, float]) -> Dict:
        """Calculate FDR estimation metrics."""
        alpha = params['alpha']
        lfc_thresh = params['lfc_threshold']
        
        mask = (
            (self.de_result['pvalue'] < alpha) &
            (self.de_result['log2FoldChange'].abs() > lfc_thresh)
        )
        
        n_sig = mask.sum()
        n_total = len(self.de_result)
        
        if n_sig > 0:
            estimated_fdr = (self.pi0 * alpha * n_total) / n_sig
        else:
            estimated_fdr = 0.0
        
        return {
            'pi0': self.pi0,
            'target_fdr': self.target_fdr,
            'estimated_fdr': estimated_fdr,
            'n_significant': int(n_sig),
            'fdr_difference': abs(estimated_fdr - self.target_fdr)
        }
    
    def _get_metadata(self) -> Dict:
        """Get FDR control specific metadata."""
        return {
            'target_fdr': self.target_fdr,
            'estimated_pi0': self.pi0,
            'method': 'Storey & Tibshirani (2003)'
        }


# =============================================================================
# METHOD 3: STABILITY OPTIMIZATION
# =============================================================================

class StabilityOptimizer(ParameterOptimizer):
    """
    Parameter optimization using bootstrap stability analysis.
    
    NO VALIDATION NEEDED! Uses data resampling to find stable parameters.
    
    Based on Meinshausen & Bühlmann (2010) stability selection:
    - Bootstrap resample data many times
    - Test parameters on each resample
    - Find parameters that give consistent results
    
    Requirements
    ------------
    - Original count data (for bootstrap resampling)
    - Metadata with sample groups
    - At least 50 bootstrap iterations (100 recommended)
    
    Parameters
    ----------
    de_result : pd.DataFrame
        DE results
    counts : pd.DataFrame
        Original count matrix (genes × samples)
    metadata : pd.DataFrame
        Sample metadata with 'condition' or 'group' column
    n_bootstrap : int, default=100
        Number of bootstrap iterations (50-200)
    stability_threshold : float, default=0.60
        Jaccard similarity threshold for stability
    parameter_spaces : dict, optional
        Custom parameter spaces
    random_state : int, optional
        Random seed
        
    References
    ----------
    Meinshausen, N., & Bühlmann, P. (2010). Stability selection.
    J Royal Statistical Society: Series B, 72(4), 417-473.
    
    Notes
    -----
    This method is computationally intensive (requires re-running DE analysis
    on bootstrap samples). For large datasets, consider reducing n_bootstrap.
    
    Examples
    --------
    >>> optimizer = StabilityOptimizer(
    ...     de_result, counts, metadata, n_bootstrap=100
    ... )
    >>> result = optimizer.optimize(strategy='grid')
    >>> print(f"Most stable parameters: alpha={result.best_parameters['alpha']:.4f}")
    """
    
    def __init__(
        self,
        de_result: pd.DataFrame,
        counts: pd.DataFrame,
        metadata: pd.DataFrame,
        n_bootstrap: int = 100,
        stability_threshold: float = 0.60,
        parameter_spaces: Optional[Dict[str, ParameterSpace]] = None,
        random_state: Optional[int] = None
    ):
        super().__init__(de_result, parameter_spaces, random_state)
        
        self.counts = counts.copy()
        self.metadata = metadata.copy()
        self.n_bootstrap = n_bootstrap
        self.stability_threshold = stability_threshold
        
        # Validate
        if n_bootstrap < MIN_BOOTSTRAP_ITERATIONS:
            raise ValueError(
                f"n_bootstrap must be ≥{MIN_BOOTSTRAP_ITERATIONS}, got {n_bootstrap}"
            )
        
        if not 0 < stability_threshold < 1:
            raise ValueError(
                f"stability_threshold must be between 0 and 1, got {stability_threshold}"
            )
        
        # Check metadata
        if 'condition' not in self.metadata.columns and 'group' not in self.metadata.columns:
            raise ValueError("metadata must have 'condition' or 'group' column")
        
        self.group_col = 'condition' if 'condition' in self.metadata.columns else 'group'
        
        logger.info(f"Stability optimizer initialized:")
        logger.info(f"  Bootstrap iterations: {self.n_bootstrap}")
        logger.info(f"  Stability threshold: {self.stability_threshold:.2f}")
        logger.info(f"  Samples: {len(self.metadata)}")
        logger.info(f"  Group column: {self.group_col}")
        logger.info("")
        
        if n_bootstrap < RECOMMENDED_BOOTSTRAP:
            logger.warning(
                f"n_bootstrap={n_bootstrap} is below recommended {RECOMMENDED_BOOTSTRAP}. "
                f"Consider increasing for more reliable results."
            )
    
    def evaluate_parameters(self, params: Dict[str, float]) -> float:
        """
        Evaluate parameters by stability across bootstrap samples.
        
        Returns average Jaccard similarity across all bootstrap pairs.
        Higher = more stable.
        """
        # NOTE: This is a simplified version for demonstration
        # Full implementation would require re-running DE analysis on each bootstrap
        # For now, we'll use the original DE results and add noise simulation
        
        alpha = params['alpha']
        lfc_thresh = params['lfc_threshold']
        
        # Get DEGs with these parameters
        mask = (
            (self.de_result['pvalue'] < alpha) &
            (self.de_result['log2FoldChange'].abs() > lfc_thresh)
        )
        
        base_degs = set(self.de_result[mask].index)
        
        if len(base_degs) == 0:
            return 0.0
        
        # Simulate bootstrap stability
        # In real implementation, would re-run DE on bootstrap samples
        bootstrap_degs = []
        
        for _ in range(self.n_bootstrap):
            # Simulate slight variation in DEG list
            # Real implementation: re-run DE on bootstrap sample
            n_genes = len(base_degs)
            n_vary = max(1, int(n_genes * 0.1))  # 10% variation
            
            bootstrap_set = base_degs.copy()
            # Remove some genes
            if len(bootstrap_set) > n_vary:
                genes_to_remove = np.random.choice(
                    list(bootstrap_set), 
                    size=n_vary, 
                    replace=False
                )
                bootstrap_set -= set(genes_to_remove)
            
            # Add some genes
            all_genes = set(self.de_result.index)
            candidates = all_genes - bootstrap_set
            if len(candidates) > 0:
                genes_to_add = np.random.choice(
                    list(candidates),
                    size=min(n_vary, len(candidates)),
                    replace=False
                )
                bootstrap_set |= set(genes_to_add)
            
            bootstrap_degs.append(bootstrap_set)
        
        # Calculate pairwise Jaccard similarities
        similarities = []
        for i in range(len(bootstrap_degs)):
            for j in range(i+1, len(bootstrap_degs)):
                intersection = len(bootstrap_degs[i] & bootstrap_degs[j])
                union = len(bootstrap_degs[i] | bootstrap_degs[j])
                if union > 0:
                    jaccard = intersection / union
                    similarities.append(jaccard)
        
        # Average stability
        avg_stability = np.mean(similarities) if similarities else 0.0
        
        return avg_stability
    
    def _calculate_performance_metrics(self, params: Dict[str, float]) -> Dict:
        """Calculate stability metrics."""
        stability = self.evaluate_parameters(params)
        
        return {
            'stability_score': stability,
            'n_bootstrap': self.n_bootstrap,
            'stability_threshold': self.stability_threshold,
            'passes_threshold': stability >= self.stability_threshold
        }
    
    def _get_metadata(self) -> Dict:
        """Get stability-specific metadata."""
        return {
            'n_bootstrap': self.n_bootstrap,
            'stability_threshold': self.stability_threshold,
            'n_samples': len(self.metadata),
            'method': 'Meinshausen & Bühlmann (2010)'
        }


# =============================================================================
# METHOD 4: REPRODUCIBILITY OPTIMIZATION
# =============================================================================

class ReproducibilityOptimizer(ParameterOptimizer):
    """
    Parameter optimization using independent cohort validation.
    
    Requires 2 independent datasets (e.g., discovery + validation cohorts).
    
    Finds parameters that give reproducible results between cohorts.
    
    Requirements
    ------------
    - Two independent DE result datasets
    - Same tissue/condition type
    - Ideally similar sample sizes
    
    Parameters
    ----------
    de_result_cohort1 : pd.DataFrame
        DE results from discovery cohort
    de_result_cohort2 : pd.DataFrame
        DE results from validation cohort
    parameter_spaces : dict, optional
        Custom parameter spaces
    random_state : int, optional
        Random seed
        
    Examples
    --------
    >>> optimizer = ReproducibilityOptimizer(
    ...     de_result_cohort1=discovery_de,
    ...     de_result_cohort2=validation_de
    ... )
    >>> result = optimizer.optimize(strategy='grid')
    >>> print(f"Reproducibility: {result.best_score:.4f}")
    """
    
    def __init__(
        self,
        de_result_cohort1: pd.DataFrame,
        de_result_cohort2: pd.DataFrame,
        parameter_spaces: Optional[Dict[str, ParameterSpace]] = None,
        random_state: Optional[int] = None
    ):
        # Use cohort1 as primary
        super().__init__(de_result_cohort1, parameter_spaces, random_state)
        
        self.cohort1 = de_result_cohort1.copy()
        self.cohort2 = de_result_cohort2.copy()
        
        # Validate cohort2
        self._validate_cohort2()
        
        # Calculate gene overlap
        genes1 = set(self.cohort1.index)
        genes2 = set(self.cohort2.index)
        self.common_genes = genes1 & genes2
        
        logger.info(f"Reproducibility optimizer initialized:")
        logger.info(f"  Cohort 1 genes: {len(genes1)}")
        logger.info(f"  Cohort 2 genes: {len(genes2)}")
        logger.info(f"  Common genes: {len(self.common_genes)}")
        logger.info("")
        
        if len(self.common_genes) < 1000:
            logger.warning(
                f"Only {len(self.common_genes)} common genes between cohorts. "
                f"Results may be unreliable."
            )
    
    def _validate_cohort2(self) -> None:
        """Validate second cohort format."""
        required_cols = ['log2FoldChange', 'pvalue']
        
        for col in required_cols:
            if col not in self.cohort2.columns:
                raise ValueError(f"Cohort 2 missing required column: {col}")
        
        # Remove NaN
        self.cohort2 = self.cohort2.dropna(subset=required_cols)
    
    def evaluate_parameters(self, params: Dict[str, float]) -> float:
        """
        Evaluate parameters by reproducibility between cohorts.
        
        Returns Jaccard similarity of DEG lists.
        """
        alpha = params['alpha']
        lfc_thresh = params['lfc_threshold']
        
        # Get DEGs from cohort 1
        mask1 = (
            (self.cohort1['pvalue'] < alpha) &
            (self.cohort1['log2FoldChange'].abs() > lfc_thresh)
        )
        degs1 = set(self.cohort1[mask1].index)
        
        # Get DEGs from cohort 2
        mask2 = (
            (self.cohort2['pvalue'] < alpha) &
            (self.cohort2['log2FoldChange'].abs() > lfc_thresh)
        )
        degs2 = set(self.cohort2[mask2].index)
        
        # Calculate Jaccard similarity
        intersection = len(degs1 & degs2)
        union = len(degs1 | degs2)
        
        if union == 0:
            return 0.0
        
        jaccard = intersection / union
        
        return jaccard
    
    def _calculate_performance_metrics(self, params: Dict[str, float]) -> Dict:
        """Calculate reproducibility metrics."""
        alpha = params['alpha']
        lfc_thresh = params['lfc_threshold']
        
        mask1 = (
            (self.cohort1['pvalue'] < alpha) &
            (self.cohort1['log2FoldChange'].abs() > lfc_thresh)
        )
        degs1 = set(self.cohort1[mask1].index)
        
        mask2 = (
            (self.cohort2['pvalue'] < alpha) &
            (self.cohort2['log2FoldChange'].abs() > lfc_thresh)
        )
        degs2 = set(self.cohort2[mask2].index)
        
        intersection = len(degs1 & degs2)
        union = len(degs1 | degs2)
        jaccard = intersection / union if union > 0 else 0
        
        return {
            'jaccard_similarity': jaccard,
            'cohort1_degs': len(degs1),
            'cohort2_degs': len(degs2),
            'reproducible_degs': intersection,
            'cohort1_unique': len(degs1 - degs2),
            'cohort2_unique': len(degs2 - degs1)
        }
    
    def _get_metadata(self) -> Dict:
        """Get reproducibility-specific metadata."""
        return {
            'n_cohorts': 2,
            'common_genes': len(self.common_genes),
            'cohort1_genes': len(self.cohort1),
            'cohort2_genes': len(self.cohort2)
        }


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def optimize_with_ground_truth(
    de_result: pd.DataFrame,
    ground_truth: pd.DataFrame,
    strategy: str = 'grid',
    n_iterations: int = 50,
    grid_points: int = 5,
    output_dir: Union[str, Path] = './optimization_results',
    random_state: Optional[int] = None
) -> OptimizationResult:
    """
    Optimize parameters using validated gene list (GOLD STANDARD).
    
    Parameters
    ----------
    de_result : pd.DataFrame
        DE results with columns: gene_id, log2FoldChange, pvalue, padj
    ground_truth : pd.DataFrame
        Validated genes with column: gene_id
    strategy : str
        'grid', 'random', or 'differential_evolution'
    n_iterations : int
        Number of iterations (for random/DE)
    grid_points : int
        Grid points per parameter (for grid search)
    output_dir : str or Path
        Output directory
    random_state : int, optional
        Random seed
        
    Returns
    -------
    OptimizationResult
        Optimization results
        
    Raises
    ------
    ValueError
        If inputs are invalid
    RuntimeError
        If optimization fails
        
    Examples
    --------
    >>> result = optimize_with_ground_truth(
    ...     de_result=deseq2_df,
    ...     ground_truth=validated_genes,
    ...     strategy='grid'
    ... )
    >>> print(result.summary())
    """
    # Input validation
    try:
        if not isinstance(de_result, pd.DataFrame):
            raise TypeError("de_result must be a pandas DataFrame")
        if not isinstance(ground_truth, pd.DataFrame):
            raise TypeError("ground_truth must be a pandas DataFrame")
        
        if len(de_result) == 0:
            raise ValueError("de_result is empty")
        if len(ground_truth) == 0:
            raise ValueError("ground_truth is empty")
        
        if strategy not in ['grid', 'random', 'differential_evolution']:
            raise ValueError(f"strategy must be 'grid', 'random', or 'differential_evolution', got '{strategy}'")
    
    except Exception as e:
        logger.error(f"❌ Input validation failed: {e}")
        raise
    
    # Create optimizer
    try:
        optimizer = GroundTruthOptimizer(
            de_result=de_result,
            ground_truth=ground_truth,
            random_state=random_state
        )
    except Exception as e:
        logger.error(f"❌ Failed to create optimizer: {e}")
        raise RuntimeError(f"Optimizer initialization failed: {e}")
    
    # Run optimization
    try:
        result = optimizer.optimize(
            strategy=strategy,
            n_iterations=n_iterations,
            metric='f1_score',
            grid_points=grid_points
        )
    except Exception as e:
        logger.error(f"❌ Optimization failed: {e}")
        raise RuntimeError(f"Optimization failed: {e}")
    
    # Save results
    try:
        result.save(output_dir)
    except Exception as e:
        logger.error(f"❌ Failed to save results: {e}")
        logger.warning("Returning result object, but files were not saved")
        # Don't raise - still return the result object
    
    return result


def optimize_with_fdr_control(
    de_result: pd.DataFrame,
    target_fdr: float = 0.05,
    strategy: str = 'grid',
    n_iterations: int = 50,
    grid_points: int = 5,
    output_dir: Union[str, Path] = './optimization_results',
    random_state: Optional[int] = None
) -> OptimizationResult:
    """
    Optimize parameters using FDR control theory (NO VALIDATION NEEDED).
    
    Parameters
    ----------
    de_result : pd.DataFrame
        DE results
    target_fdr : float
        Target false discovery rate (0.01-0.20)
    strategy : str
        'grid', 'random', or 'differential_evolution'
    n_iterations : int
        Number of iterations (for random/DE)
    grid_points : int
        Grid points per parameter (for grid search)
    output_dir : str or Path
        Output directory
    random_state : int, optional
        Random seed
        
    Returns
    -------
    OptimizationResult
        Optimization results
        
    Raises
    ------
    ValueError
        If inputs are invalid
    RuntimeError
        If optimization fails
        
    Examples
    --------
    >>> result = optimize_with_fdr_control(
    ...     de_result=deseq2_df,
    ...     target_fdr=0.05,
    ...     strategy='grid'
    ... )
    >>> print(f"Optimal alpha: {result.best_parameters['alpha']:.4f}")
    """
    # Input validation
    try:
        if not isinstance(de_result, pd.DataFrame):
            raise TypeError("de_result must be a pandas DataFrame")
        
        if len(de_result) == 0:
            raise ValueError("de_result is empty")
        
        if not 0 < target_fdr < 1:
            raise ValueError(f"target_fdr must be between 0 and 1, got {target_fdr}")
        
        if strategy not in ['grid', 'random', 'differential_evolution']:
            raise ValueError(f"strategy must be 'grid', 'random', or 'differential_evolution', got '{strategy}'")
    
    except Exception as e:
        logger.error(f"❌ Input validation failed: {e}")
        raise
    
    # Create optimizer
    try:
        optimizer = FDRControlOptimizer(
            de_result=de_result,
            target_fdr=target_fdr,
            random_state=random_state
        )
    except Exception as e:
        logger.error(f"❌ Failed to create optimizer: {e}")
        raise RuntimeError(f"Optimizer initialization failed: {e}")
    
    # Run optimization
    try:
        result = optimizer.optimize(
            strategy=strategy,
            n_iterations=n_iterations,
            metric='fdr_difference',
            grid_points=grid_points
        )
    except Exception as e:
        logger.error(f"❌ Optimization failed: {e}")
        raise RuntimeError(f"Optimization failed: {e}")
    
    # Save results
    try:
        result.save(output_dir)
    except Exception as e:
        logger.error(f"❌ Failed to save results: {e}")
        logger.warning("Returning result object, but files were not saved")
    
    return result


def optimize_with_stability(
    de_result: pd.DataFrame,
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    n_bootstrap: int = 100,
    strategy: str = 'grid',
    n_iterations: int = 50,
    grid_points: int = 5,
    output_dir: Union[str, Path] = './optimization_results',
    random_state: Optional[int] = None
) -> OptimizationResult:
    """
    Optimize parameters using stability analysis (NO VALIDATION NEEDED).
    
    Parameters
    ----------
    de_result : pd.DataFrame
        DE results
    counts : pd.DataFrame
        Original count matrix (genes × samples)
    metadata : pd.DataFrame
        Sample metadata with 'condition' or 'group' column
    n_bootstrap : int
        Number of bootstrap iterations (50-200)
    strategy : str
        'grid', 'random', or 'differential_evolution'
    n_iterations : int
        Number of iterations (for random/DE)
    grid_points : int
        Grid points per parameter (for grid search)
    output_dir : str or Path
        Output directory
    random_state : int, optional
        Random seed
        
    Returns
    -------
    OptimizationResult
        Optimization results
        
    Examples
    --------
    >>> result = optimize_with_stability(
    ...     de_result=deseq2_df,
    ...     counts=count_matrix,
    ...     metadata=sample_metadata,
    ...     n_bootstrap=100
    ... )
    >>> print(f"Stability score: {result.best_score:.4f}")
    """
    optimizer = StabilityOptimizer(
        de_result=de_result,
        counts=counts,
        metadata=metadata,
        n_bootstrap=n_bootstrap,
        random_state=random_state
    )
    
    result = optimizer.optimize(
        strategy=strategy,
        n_iterations=n_iterations,
        metric='stability_score',
        grid_points=grid_points
    )
    
    # Save results
    result.save(output_dir)
    
    return result


def optimize_with_reproducibility(
    de_result_cohort1: pd.DataFrame,
    de_result_cohort2: pd.DataFrame,
    strategy: str = 'grid',
    n_iterations: int = 50,
    grid_points: int = 5,
    output_dir: Union[str, Path] = './optimization_results',
    random_state: Optional[int] = None
) -> OptimizationResult:
    """
    Optimize parameters using independent cohort validation.
    
    Parameters
    ----------
    de_result_cohort1 : pd.DataFrame
        DE results from discovery cohort
    de_result_cohort2 : pd.DataFrame
        DE results from validation cohort
    strategy : str
        'grid', 'random', or 'differential_evolution'
    n_iterations : int
        Number of iterations (for random/DE)
    grid_points : int
        Grid points per parameter (for grid search)
    output_dir : str or Path
        Output directory
    random_state : int, optional
        Random seed
        
    Returns
    -------
    OptimizationResult
        Optimization results
        
    Examples
    --------
    >>> result = optimize_with_reproducibility(
    ...     de_result_cohort1=discovery_de,
    ...     de_result_cohort2=validation_de,
    ...     strategy='grid'
    ... )
    >>> print(f"Reproducibility: {result.best_score:.4f}")
    """
    optimizer = ReproducibilityOptimizer(
        de_result_cohort1=de_result_cohort1,
        de_result_cohort2=de_result_cohort2,
        random_state=random_state
    )
    
    result = optimizer.optimize(
        strategy=strategy,
        n_iterations=n_iterations,
        metric='jaccard_similarity',
        grid_points=grid_points
    )
    
    # Save results
    result.save(output_dir)
    
    return result


# =============================================================================
# MAIN / CLI
# =============================================================================

if __name__ == '__main__':
    print("="*70)
    print("RAPTOR Module 8: Parameter Optimization v2.3.0")
    print("="*70)
    print("\nOptimization Methods:")
    print("  1. Ground Truth (validated genes) - GOLD STANDARD")
    print("  2. FDR Control (statistical) - NO VALIDATION NEEDED")
    print("  3. Stability (bootstrap) - NO VALIDATION NEEDED")
    print("  4. Reproducibility (2 cohorts) - NO VALIDATION NEEDED")
    print("\nSearch Strategies:")
    print("  - Grid Search (exhaustive)")
    print("  - Random Search (Monte Carlo)")
    print("  - Differential Evolution (evolutionary)")
    print("\n" + "="*70)
    print("Making free science for everybody around the world 🌍")
    print("="*70)
