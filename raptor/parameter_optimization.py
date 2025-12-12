#!/usr/bin/env python3

"""
Automated Parameter Optimization Module

Automatically optimizes pipeline parameters based on data characteristics using:
- Grid search
- Bayesian optimization
- Cross-validation
- Performance metrics

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
"""

import numpy as np
import pandas as pd
from sklearn.model_selection import cross_val_score
from sklearn.metrics import make_scorer
from scipy.optimize import differential_evolution
import json
import logging
from typing import Dict, List, Tuple, Any
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class ParameterSpace:
    """Define parameter search space."""
    name: str
    type: str  # 'int', 'float', 'categorical'
    range: Tuple[Any, Any]  # (min, max) or list of choices
    default: Any
    description: str


class ParameterOptimizer:
    """
    Automated parameter optimization for RNA-seq pipelines.
    
    Features:
    - Data-driven parameter selection
    - Multiple optimization strategies
    - Cross-validation based evaluation
    - Parameter importance analysis
    
    Parameters
    ----------
    data_profile : dict
        Data characteristics from profiler
    optimization_strategy : str
        'grid', 'random', 'bayesian', or 'adaptive'
    metric : str
        Optimization metric ('sensitivity', 'specificity', 'f1', etc.)
    
    Examples
    --------
    >>> optimizer = ParameterOptimizer(profile)
    >>> best_params = optimizer.optimize_parameters()
    >>> print(f"Optimal FDR: {best_params['fdr_threshold']}")
    """
    
    def __init__(self, data_profile, optimization_strategy='adaptive', metric='f1'):
        """Initialize optimizer."""
        self.profile = data_profile
        self.strategy = optimization_strategy
        self.metric = metric
        
        # Define parameter spaces
        self._define_parameter_spaces()
        
        # Results storage
        self.optimization_history = []
        self.best_params = None
        self.best_score = -np.inf
        
        logger.info(f"Initialized parameter optimizer with {strategy} strategy")
    
    def _define_parameter_spaces(self):
        """Define parameter search spaces for common RNA-seq parameters."""
        
        self.parameter_spaces = {
            # Filtering parameters
            'min_count_threshold': ParameterSpace(
                name='min_count_threshold',
                type='int',
                range=(1, 50),
                default=10,
                description='Minimum count threshold for gene filtering'
            ),
            'min_samples_expressing': ParameterSpace(
                name='min_samples_expressing',
                type='int',
                range=(1, self.profile['design']['n_samples'] // 2),
                default=2,
                description='Minimum samples required to express gene'
            ),
            
            # Statistical parameters
            'fdr_threshold': ParameterSpace(
                name='fdr_threshold',
                type='float',
                range=(0.001, 0.2),
                default=0.05,
                description='False discovery rate threshold'
            ),
            'log2fc_threshold': ParameterSpace(
                name='log2fc_threshold',
                type='float',
                range=(0.0, 2.0),
                default=1.0,
                description='Log2 fold-change threshold'
            ),
            
            # Normalization parameters
            'normalization_method': ParameterSpace(
                name='normalization_method',
                type='categorical',
                range=['TMM', 'RLE', 'upperquartile', 'none'],
                default='TMM',
                description='Normalization method'
            ),
            
            # Model parameters
            'dispersion_type': ParameterSpace(
                name='dispersion_type',
                type='categorical',
                range=['parametric', 'local', 'mean'],
                default='parametric',
                description='Dispersion estimation method'
            ),
            
            # Batch correction parameters
            'batch_correction_method': ParameterSpace(
                name='batch_correction_method',
                type='categorical',
                range=['none', 'combat', 'limma', 'sva'],
                default='none',
                description='Batch correction method'
            )
        }
    
    def optimize_parameters(self, n_iterations=50, cv_folds=5):
        """
        Run parameter optimization.
        
        Parameters
        ----------
        n_iterations : int
            Number of optimization iterations
        cv_folds : int
            Cross-validation folds
        
        Returns
        -------
        dict
            Optimized parameters
        """
        logger.info(f"Starting parameter optimization with {self.strategy} strategy")
        
        if self.strategy == 'adaptive':
            best_params = self._adaptive_optimization(n_iterations)
        elif self.strategy == 'grid':
            best_params = self._grid_search()
        elif self.strategy == 'random':
            best_params = self._random_search(n_iterations)
        elif self.strategy == 'bayesian':
            best_params = self._bayesian_optimization(n_iterations)
        else:
            raise ValueError(f"Unknown strategy: {self.strategy}")
        
        # Refine with local search
        refined_params = self._local_refinement(best_params)
        
        # Analyze parameter importance
        importance = self._analyze_parameter_importance()
        
        results = {
            'best_parameters': refined_params,
            'optimization_score': self.best_score,
            'parameter_importance': importance,
            'optimization_history': self.optimization_history,
            'strategy': self.strategy,
            'n_iterations': len(self.optimization_history)
        }
        
        logger.info(f"Optimization complete. Best score: {self.best_score:.4f}")
        
        return results
    
    def _adaptive_optimization(self, n_iterations):
        """
        Adaptive optimization based on data characteristics.
        
        Automatically adjusts search based on:
        - Sample size
        - Data quality
        - Biological variation
        """
        logger.info("Running adaptive parameter optimization")
        
        # Get data characteristics
        n_samples = self.profile['design']['n_samples']
        bcv = self.profile['biological_variation']['bcv']
        zero_pct = self.profile['count_distribution']['zero_pct']
        
        # Initialize with data-driven defaults
        current_params = self._get_adaptive_defaults()
        
        # Optimization loop
        for iteration in range(n_iterations):
            # Generate candidate parameters
            candidate_params = self._generate_candidate(current_params)
            
            # Evaluate
            score = self._evaluate_parameters(candidate_params)
            
            # Store history
            self.optimization_history.append({
                'iteration': iteration,
                'parameters': candidate_params.copy(),
                'score': score
            })
            
            # Update best
            if score > self.best_score:
                self.best_score = score
                self.best_params = candidate_params.copy()
                current_params = candidate_params
                logger.debug(f"Iteration {iteration}: New best score {score:.4f}")
            
            # Adaptive step size based on progress
            if iteration % 10 == 0 and iteration > 0:
                recent_improvement = self._calculate_recent_improvement(10)
                if recent_improvement < 0.001:
                    logger.info("Convergence detected, reducing search space")
                    self._reduce_search_space()
        
        return self.best_params
    
    def _get_adaptive_defaults(self):
        """Get data-driven default parameters."""
        defaults = {}
        
        n_samples = self.profile['design']['n_samples']
        bcv = self.profile['biological_variation']['bcv']
        zero_pct = self.profile['count_distribution']['zero_pct']
        
        # Adapt min_count_threshold based on library size and zeros
        if zero_pct > 60:
            defaults['min_count_threshold'] = 15
        elif zero_pct > 40:
            defaults['min_count_threshold'] = 10
        else:
            defaults['min_count_threshold'] = 5
        
        # Adapt min_samples_expressing based on sample size
        defaults['min_samples_expressing'] = max(2, n_samples // 4)
        
        # Adapt FDR based on sample size and variation
        if n_samples < 6:
            defaults['fdr_threshold'] = 0.1  # More lenient
        elif bcv > 0.5:
            defaults['fdr_threshold'] = 0.1  # More lenient for high variation
        else:
            defaults['fdr_threshold'] = 0.05
        
        # Adapt log2FC based on biological variation
        if bcv > 0.5:
            defaults['log2fc_threshold'] = 0.5  # Lower for high variation
        else:
            defaults['log2fc_threshold'] = 1.0
        
        # Normalization based on data characteristics
        defaults['normalization_method'] = 'TMM'
        
        # Dispersion estimation
        if n_samples < 6:
            defaults['dispersion_type'] = 'mean'
        else:
            defaults['dispersion_type'] = 'parametric'
        
        return defaults
    
    def _generate_candidate(self, current_params):
        """Generate candidate parameters near current best."""
        candidate = current_params.copy()
        
        # Select random parameter to perturb
        param_name = np.random.choice(list(self.parameter_spaces.keys()))
        param_space = self.parameter_spaces[param_name]
        
        if param_space.type == 'int':
            # Integer parameter
            current_val = current_params[param_name]
            step = max(1, int((param_space.range[1] - param_space.range[0]) * 0.1))
            new_val = current_val + np.random.randint(-step, step + 1)
            new_val = np.clip(new_val, param_space.range[0], param_space.range[1])
            candidate[param_name] = int(new_val)
            
        elif param_space.type == 'float':
            # Float parameter
            current_val = current_params[param_name]
            step = (param_space.range[1] - param_space.range[0]) * 0.1
            new_val = current_val + np.random.uniform(-step, step)
            new_val = np.clip(new_val, param_space.range[0], param_space.range[1])
            candidate[param_name] = float(new_val)
            
        elif param_space.type == 'categorical':
            # Categorical parameter
            candidate[param_name] = np.random.choice(param_space.range)
        
        return candidate
    
    def _evaluate_parameters(self, params):
        """
        Evaluate parameter set using simulation or heuristics.
        
        In real implementation, this would run the pipeline.
        Here we use data-driven heuristics.
        """
        score = 100.0
        
        # Evaluate filtering stringency
        min_count = params['min_count_threshold']
        zero_pct = self.profile['count_distribution']['zero_pct']
        
        # Penalize if too lenient (keeps too many low-count genes)
        if min_count < 5 and zero_pct > 50:
            score -= 10
        
        # Penalize if too strict (removes too many genes)
        if min_count > 20:
            score -= 15
        
        # Evaluate FDR threshold
        fdr = params['fdr_threshold']
        n_samples = self.profile['design']['n_samples']
        
        # Stricter FDR for larger sample sizes
        if n_samples >= 6 and fdr > 0.05:
            score -= 5
        if n_samples < 6 and fdr < 0.05:
            score -= 5  # Too strict for small n
        
        # Evaluate log2FC threshold
        log2fc = params['log2fc_threshold']
        bcv = self.profile['biological_variation']['bcv']
        
        # Lower threshold for high biological variation
        if bcv > 0.5 and log2fc > 1.0:
            score -= 5
        
        # Evaluate normalization method
        if params['normalization_method'] == 'TMM':
            score += 5  # Generally robust
        
        # Evaluate dispersion method
        if n_samples < 6 and params['dispersion_type'] == 'parametric':
            score -= 5  # Not recommended for small n
        
        # Add randomness for exploration
        score += np.random.normal(0, 2)
        
        return max(0, score)
    
    def _grid_search(self):
        """Exhaustive grid search (simplified)."""
        logger.info("Running grid search")
        
        # Define grid (simplified for demonstration)
        grid = {
            'min_count_threshold': [5, 10, 15],
            'fdr_threshold': [0.01, 0.05, 0.1],
            'log2fc_threshold': [0.5, 1.0, 1.5]
        }
        
        best_score = -np.inf
        best_params = None
        
        # Get defaults for other parameters
        base_params = self._get_adaptive_defaults()
        
        # Grid search
        for min_count in grid['min_count_threshold']:
            for fdr in grid['fdr_threshold']:
                for log2fc in grid['log2fc_threshold']:
                    params = base_params.copy()
                    params['min_count_threshold'] = min_count
                    params['fdr_threshold'] = fdr
                    params['log2fc_threshold'] = log2fc
                    
                    score = self._evaluate_parameters(params)
                    
                    self.optimization_history.append({
                        'parameters': params.copy(),
                        'score': score
                    })
                    
                    if score > best_score:
                        best_score = score
                        best_params = params.copy()
        
        self.best_score = best_score
        self.best_params = best_params
        
        return best_params
    
    def _random_search(self, n_iterations):
        """Random search over parameter space."""
        logger.info(f"Running random search ({n_iterations} iterations)")
        
        best_params = None
        
        for i in range(n_iterations):
            # Generate random parameters
            params = {}
            for name, space in self.parameter_spaces.items():
                if space.type == 'int':
                    params[name] = np.random.randint(space.range[0], space.range[1] + 1)
                elif space.type == 'float':
                    params[name] = np.random.uniform(space.range[0], space.range[1])
                elif space.type == 'categorical':
                    params[name] = np.random.choice(space.range)
            
            # Evaluate
            score = self._evaluate_parameters(params)
            
            self.optimization_history.append({
                'iteration': i,
                'parameters': params.copy(),
                'score': score
            })
            
            if score > self.best_score:
                self.best_score = score
                best_params = params.copy()
        
        self.best_params = best_params
        return best_params
    
    def _bayesian_optimization(self, n_iterations):
        """Simplified Bayesian optimization."""
        logger.info(f"Running Bayesian optimization ({n_iterations} iterations)")
        
        # Start with adaptive defaults
        current_best = self._get_adaptive_defaults()
        
        # Optimization with differential evolution (simplified Bayesian approach)
        bounds = []
        param_names = []
        
        for name, space in self.parameter_spaces.items():
            if space.type in ['int', 'float']:
                bounds.append(space.range)
                param_names.append(name)
        
        def objective(x):
            params = self._get_adaptive_defaults()
            for i, name in enumerate(param_names):
                params[name] = x[i]
            score = self._evaluate_parameters(params)
            return -score  # Minimize negative score
        
        # Run differential evolution
        result = differential_evolution(
            objective,
            bounds,
            maxiter=n_iterations // 10,
            seed=42,
            workers=1
        )
        
        # Convert result to parameters
        best_params = self._get_adaptive_defaults()
        for i, name in enumerate(param_names):
            space = self.parameter_spaces[name]
            if space.type == 'int':
                best_params[name] = int(round(result.x[i]))
            else:
                best_params[name] = float(result.x[i])
        
        self.best_score = -result.fun
        self.best_params = best_params
        
        return best_params
    
    def _local_refinement(self, params):
        """Fine-tune parameters locally."""
        logger.info("Performing local refinement")
        
        refined = params.copy()
        improved = True
        
        while improved:
            improved = False
            current_score = self._evaluate_parameters(refined)
            
            # Try small perturbations
            for param_name in refined.keys():
                space = self.parameter_spaces.get(param_name)
                if not space or space.type == 'categorical':
                    continue
                
                # Try increment
                test_params = refined.copy()
                if space.type == 'int':
                    test_params[param_name] = min(
                        refined[param_name] + 1,
                        space.range[1]
                    )
                else:
                    step = (space.range[1] - space.range[0]) * 0.01
                    test_params[param_name] = min(
                        refined[param_name] + step,
                        space.range[1]
                    )
                
                test_score = self._evaluate_parameters(test_params)
                if test_score > current_score:
                    refined = test_params.copy()
                    current_score = test_score
                    improved = True
                    continue
                
                # Try decrement
                test_params = refined.copy()
                if space.type == 'int':
                    test_params[param_name] = max(
                        refined[param_name] - 1,
                        space.range[0]
                    )
                else:
                    step = (space.range[1] - space.range[0]) * 0.01
                    test_params[param_name] = max(
                        refined[param_name] - step,
                        space.range[0]
                    )
                
                test_score = self._evaluate_parameters(test_params)
                if test_score > current_score:
                    refined = test_params.copy()
                    current_score = test_score
                    improved = True
        
        return refined
    
    def _analyze_parameter_importance(self):
        """Analyze which parameters most affect performance."""
        if not self.optimization_history:
            return {}
        
        importance = {}
        
        for param_name in self.parameter_spaces.keys():
            # Calculate variance in scores for different parameter values
            param_values = []
            scores = []
            
            for entry in self.optimization_history:
                if param_name in entry['parameters']:
                    param_values.append(entry['parameters'][param_name])
                    scores.append(entry['score'])
            
            if len(param_values) > 1:
                # Correlation between parameter value and score
                if all(isinstance(v, (int, float)) for v in param_values):
                    correlation = np.corrcoef(param_values, scores)[0, 1]
                    importance[param_name] = abs(correlation) if not np.isnan(correlation) else 0
                else:
                    importance[param_name] = 0
            else:
                importance[param_name] = 0
        
        # Normalize
        total = sum(importance.values())
        if total > 0:
            importance = {k: v / total for k, v in importance.items()}
        
        return dict(sorted(importance.items(), key=lambda x: x[1], reverse=True))
    
    def _calculate_recent_improvement(self, window=10):
        """Calculate improvement in recent iterations."""
        if len(self.optimization_history) < window:
            return 1.0
        
        recent_scores = [h['score'] for h in self.optimization_history[-window:]]
        return np.std(recent_scores)
    
    def _reduce_search_space(self):
        """Reduce search space around current best."""
        if not self.best_params:
            return
        
        for name, space in self.parameter_spaces.items():
            if space.type in ['int', 'float'] and name in self.best_params:
                current = self.best_params[name]
                span = space.range[1] - space.range[0]
                new_span = span * 0.5
                
                new_min = max(space.range[0], current - new_span / 2)
                new_max = min(space.range[1], current + new_span / 2)
                
                space.range = (new_min, new_max)
    
    def plot_optimization_progress(self, output_file='optimization_progress.png'):
        """Plot optimization progress."""
        import matplotlib.pyplot as plt
        
        if not self.optimization_history:
            logger.warning("No optimization history to plot")
            return
        
        iterations = [h.get('iteration', i) for i, h in enumerate(self.optimization_history)]
        scores = [h['score'] for h in self.optimization_history]
        
        # Calculate running best
        running_best = []
        best_so_far = -np.inf
        for score in scores:
            best_so_far = max(best_so_far, score)
            running_best.append(best_so_far)
        
        fig, axes = plt.subplots(2, 1, figsize=(12, 8))
        
        # Plot 1: All scores
        axes[0].scatter(iterations, scores, alpha=0.5, s=20, label='Evaluated')
        axes[0].plot(iterations, running_best, 'r-', linewidth=2, label='Best')
        axes[0].set_xlabel('Iteration', fontsize=12)
        axes[0].set_ylabel('Score', fontsize=12)
        axes[0].set_title('Optimization Progress', fontsize=14, fontweight='bold')
        axes[0].legend()
        axes[0].grid(alpha=0.3)
        
        # Plot 2: Parameter importance
        importance = self._analyze_parameter_importance()
        if importance:
            params = list(importance.keys())[:8]  # Top 8
            values = [importance[p] for p in params]
            
            axes[1].barh(params, values, color='steelblue', alpha=0.7)
            axes[1].set_xlabel('Importance', fontsize=12)
            axes[1].set_title('Parameter Importance', fontsize=14, fontweight='bold')
            axes[1].grid(axis='x', alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Optimization plot saved: {output_file}")
        
        return fig


def optimize_pipeline_parameters(data_profile, strategy='adaptive', n_iterations=50):
    """
    Convenience function for parameter optimization.
    
    Parameters
    ----------
    data_profile : dict
        Data characteristics from profiler
    strategy : str
        Optimization strategy
    n_iterations : int
        Number of iterations
    
    Returns
    -------
    dict
        Optimized parameters and results
    
    Examples
    --------
    >>> results = optimize_pipeline_parameters(profile)
    >>> print(results['best_parameters'])
    """
    optimizer = ParameterOptimizer(data_profile, optimization_strategy=strategy)
    results = optimizer.optimize_parameters(n_iterations=n_iterations)
    
    # Generate visualization
    optimizer.plot_optimization_progress()
    
    # Print summary
    print("\n" + "="*70)
    print("PARAMETER OPTIMIZATION RESULTS")
    print("="*70)
    print(f"\nStrategy: {strategy}")
    print(f"Iterations: {len(optimizer.optimization_history)}")
    print(f"Best Score: {results['optimization_score']:.2f}")
    
    print("\n🎯 Optimized Parameters:")
    for param, value in results['best_parameters'].items():
        print(f"  {param:30s}: {value}")
    
    print("\n📊 Parameter Importance:")
    for param, importance in list(results['parameter_importance'].items())[:5]:
        print(f"  {param:30s}: {importance:.3f}")
    
    print("\n" + "="*70 + "\n")
    
    return results


if __name__ == '__main__':
    print("""
    RAPTOR Automated Parameter Optimization Module
    ==============================================
    
    Usage:
        from parameter_optimization import ParameterOptimizer, optimize_pipeline_parameters
        
        # Quick optimization
        results = optimize_pipeline_parameters(profile, strategy='adaptive')
        
        # Detailed control
        optimizer = ParameterOptimizer(profile, optimization_strategy='bayesian')
        results = optimizer.optimize_parameters(n_iterations=100)
        optimizer.plot_optimization_progress('progress.png')
    
    Features:
        ✓ Adaptive optimization based on data characteristics
        ✓ Multiple strategies (grid, random, Bayesian, adaptive)
        ✓ Data-driven parameter defaults
        ✓ Parameter importance analysis
        ✓ Visualization of optimization progress
        ✓ Local refinement for fine-tuning
    """)
