"""
ML-Based Pipeline Recommender (TRUE Machine Learning) - VALIDATED VERSION

Machine learning recommender trained on ACTUAL pipeline performance data.
Unlike heuristic approaches, this learns from real benchmark results.

Training Pipeline:
1. Generate diverse simulated datasets with known ground truth
2. Run ALL DE pipelines on each dataset
3. Evaluate against ground truth (TPR, FDR, AUC)
4. Label each dataset with the BEST performing pipeline
5. Train classifier to predict best pipeline from data features

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0 (VALIDATED)
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass, field
import logging
import json
import pickle
from datetime import datetime
import warnings

# ML libraries
try:
    from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
    from sklearn.preprocessing import StandardScaler, LabelEncoder
    from sklearn.model_selection import cross_val_score, train_test_split, StratifiedKFold
    from sklearn.metrics import classification_report, confusion_matrix, accuracy_score
    import joblib
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False
    warnings.warn("scikit-learn not available. Install with: pip install scikit-learn")

# RAPTOR v2.2.0 imports
from raptor.utils.validation import (
    validate_count_matrix,
    validate_metadata,
    validate_file_path,
    validate_directory_path,
    validate_numeric_range,
    validate_positive_integer
)

from raptor.utils.errors import (
    handle_errors,
    ValidationError,
    DependencyError,
    check_file_exists,
    validate_output_writable
)

logger = logging.getLogger(__name__)


# =============================================================================
# CLI PARAMETERS & MODULE METADATA - v2.2.0
# =============================================================================

ML_RECOMMENDER_CLI_PARAMS = {
    'profile': {
        'flag': '--profile',
        'short': '-p',
        'required': True,
        'type': 'file',
        'help': 'Data profile JSON file from profiler module',
        'validation': 'Must be valid profile JSON'
    },
    'model_file': {
        'flag': '--model',
        'short': '-m',
        'type': 'file',
        'help': 'Pre-trained ML model file (optional, will train if not provided)',
        'validation': 'Must be pickle file from MLPipelineRecommender'
    },
    'n_datasets': {
        'flag': '--n-datasets',
        'short': '-n',
        'default': 100,
        'type': 'int',
        'range': [10, 10000],
        'help': 'Number of simulated datasets for training (default: 100)'
    },
    'model_type': {
        'flag': '--model-type',
        'default': 'random_forest',
        'choices': ['random_forest', 'gradient_boosting'],
        'help': 'ML model type (default: random_forest)'
    },
    'criterion': {
        'flag': '--criterion',
        'short': '-c',
        'default': 'f1_score',
        'choices': ['f1_score', 'auc', 'tpr', 'fdr_control'],
        'help': 'Optimization criterion (default: f1_score)'
    },
    'output': {
        'flag': '--output',
        'short': '-o',
        'default': 'ml_recommendation.json',
        'type': 'file',
        'help': 'Output recommendation JSON file'
    },
    'save_model': {
        'flag': '--save-model',
        'type': 'file',
        'help': 'Save trained model to file (optional)'
    },
    'verbose': {
        'flag': '--verbose',
        'short': '-v',
        'action': 'store_true',
        'help': 'Verbose training/prediction output'
    }
}

# Module metadata  
ML_MODULE_NAME = "ML Pipeline Recommender"
ML_MODULE_VERSION = "2.2.0"
ML_MODULE_ID = "M4-ML"

__all__ = [
    # Main class
    'MLPipelineRecommender',
    
    # Data classes
    'MLRecommendation',
    'TrainingConfig',
    
    # Convenience functions
    'train_ml_recommender',
    'predict_with_ml',
    
    # CLI integration
    'ML_RECOMMENDER_CLI_PARAMS',
    'ML_MODULE_NAME',
    'ML_MODULE_VERSION',
    'ML_MODULE_ID'
]


@dataclass
class TrainingConfig:
    """Configuration for ML training.
    
    Parameters
    ----------
    n_datasets : int
        Number of synthetic datasets to generate for training (must be >= 10)
    model_type : str
        'random_forest' or 'gradient_boosting'
    criterion : str
        Performance criterion: 'f1_score', 'auc', 'tpr', 'fdr_control'
    test_size : float
        Fraction of data for testing (0-1)
    cv_folds : int
        Cross-validation folds (must be >= 2)
    n_jobs : int
        Parallel jobs (-1 for all cores)
    seed : int
        Random seed
    """
    n_datasets: int = 100
    model_type: str = 'random_forest'
    criterion: str = 'f1_score'
    test_size: float = 0.2
    cv_folds: int = 5
    n_jobs: int = -1
    seed: int = 42
    
    def __post_init__(self):
        """Validate configuration after initialization."""
        # Validate n_datasets
        if self.n_datasets < 10:
            raise ValidationError(
                f"n_datasets must be >= 10, got {self.n_datasets}. "
                "Need sufficient datasets for reliable training."
            )
        
        # Validate model_type
        valid_models = ['random_forest', 'gradient_boosting']
        if self.model_type not in valid_models:
            raise ValidationError(
                f"model_type must be one of {valid_models}, got '{self.model_type}'"
            )
        
        # Validate criterion
        valid_criteria = ['f1_score', 'auc', 'tpr', 'fdr_control']
        if self.criterion not in valid_criteria:
            raise ValidationError(
                f"criterion must be one of {valid_criteria}, got '{self.criterion}'"
            )
        
        # Validate test_size
        if not 0 < self.test_size < 1:
            raise ValidationError(
                f"test_size must be between 0 and 1, got {self.test_size}"
            )
        
        # Validate cv_folds
        if self.cv_folds < 2:
            raise ValidationError(
                f"cv_folds must be >= 2, got {self.cv_folds}"
            )


@dataclass
class MLRecommendation:
    """ML-based pipeline recommendation.
    
    Attributes
    ----------
    pipeline : str
        Recommended pipeline name
    confidence : float
        Prediction confidence (probability)
    alternatives : List[Tuple[str, float]]
        Alternative pipelines with probabilities
    feature_importances : Dict[str, float]
        Top contributing features
    model_version : str
        Version of trained model
    """
    pipeline: str
    confidence: float
    alternatives: List[Tuple[str, float]]
    feature_importances: Dict[str, float]
    model_version: str
    
    def __post_init__(self):
        """Validate recommendation after initialization."""
        # Validate confidence
        if not 0 <= self.confidence <= 1:
            raise ValidationError(
                f"confidence must be between 0 and 1, got {self.confidence}"
            )
        
        # Validate alternatives
        for alt_name, alt_prob in self.alternatives:
            if not isinstance(alt_name, str):
                raise ValidationError(
                    f"Alternative pipeline name must be string, got {type(alt_name)}"
                )
            if not 0 <= alt_prob <= 1:
                raise ValidationError(
                    f"Alternative probability must be between 0 and 1, got {alt_prob}"
                )
    
    def summary(self) -> str:
        """Generate summary string."""
        lines = [
            f"🤖 ML Recommendation: {self.pipeline}",
            f"   Confidence: {self.confidence:.1%}",
            "",
            "   Alternatives:"
        ]
        for alt, prob in self.alternatives[:3]:
            lines.append(f"     • {alt}: {prob:.1%}")
        
        lines.extend([
            "",
            "   Key factors:"
        ])
        for feat, imp in list(self.feature_importances.items())[:5]:
            lines.append(f"     • {feat}: {imp:.3f}")
        
        return "\n".join(lines)


class MLPipelineRecommender:
    """
    TRUE Machine Learning Pipeline Recommender.
    
    Trains on actual benchmark performance, not heuristics.
    
    Training Process:
    1. Generate N diverse synthetic RNA-seq datasets
    2. Run DESeq2, edgeR, limma-voom, Wilcoxon on each
    3. Evaluate each against ground truth DE genes
    4. Label dataset with best-performing pipeline
    5. Train classifier (RF or GB) on data features
    
    Parameters
    ----------
    model_type : str
        'random_forest' or 'gradient_boosting'
    model_path : str, optional
        Path to load pre-trained model
    
    Examples
    --------
    >>> # Train new model
    >>> recommender = MLPipelineRecommender()
    >>> recommender.train(n_datasets=100)
    >>> recommender.save('models/ml_recommender.pkl')
    
    >>> # Use trained model
    >>> recommender = MLPipelineRecommender(model_path='models/')
    >>> profile = profile_data_quick(counts, metadata)
    >>> rec = recommender.recommend(profile)
    >>> print(rec.pipeline, rec.confidence)
    """
    
    VERSION = "2.2.0"
    SUPPORTED_PIPELINES = ['DESeq2', 'edgeR', 'limma_voom', 'Wilcoxon']
    
    @handle_errors(exit_on_error=False)
    def __init__(
        self,
        model_type: str = 'random_forest',
        model_path: Optional[str] = None
    ):
        """Initialize ML recommender with validation."""
        # Check sklearn availability
        if not SKLEARN_AVAILABLE:
            raise DependencyError(
                "scikit-learn required for ML recommender. "
                "Install with: pip install scikit-learn joblib"
            )
        
        # Validate model_type
        valid_models = ['random_forest', 'gradient_boosting']
        if model_type not in valid_models:
            raise ValidationError(
                f"model_type must be one of {valid_models}, got '{model_type}'"
            )
        
        self.model_type = model_type
        self.model = None
        self.scaler = StandardScaler()
        self.label_encoder = LabelEncoder()
        self.feature_names = None
        self.is_trained = False
        self.training_metadata = {}
        
        if model_path:
            # Validate path exists
            path = Path(model_path)
            if not path.exists():
                raise FileNotFoundError(
                    f"Model path not found: {model_path}\n"
                    f"    Current directory: {Path.cwd()}\n"
                    f"    Train a new model with .train() or provide correct path"
                )
            self.load(model_path)
        else:
            self._initialize_model()
        
        logger.info(f"✓ Initialized ML recommender: {model_type}")
    
    def _initialize_model(self):
        """Initialize the ML model."""
        if self.model_type == 'random_forest':
            self.model = RandomForestClassifier(
                n_estimators=200,
                max_depth=15,
                min_samples_split=5,
                min_samples_leaf=2,
                max_features='sqrt',
                class_weight='balanced',
                random_state=42,
                n_jobs=-1
            )
        elif self.model_type == 'gradient_boosting':
            self.model = GradientBoostingClassifier(
                n_estimators=200,
                max_depth=5,
                learning_rate=0.1,
                min_samples_split=5,
                min_samples_leaf=2,
                random_state=42
            )
        else:
            raise ValidationError(f"Unknown model type: {self.model_type}")
    
    @handle_errors(exit_on_error=False)
    def train(
        self,
        n_datasets: int = 100,
        criterion: str = 'f1_score',
        verbose: bool = True
    ) -> Dict:
        """
        Train the ML model on synthetic benchmark data.
        
        This is the REAL ML training - actually runs pipelines and evaluates!
        
        Parameters
        ----------
        n_datasets : int
            Number of synthetic datasets to generate (must be >= 10)
        criterion : str
            How to determine "best" pipeline: 'f1_score', 'auc', 'tpr', 'fdr_control'
        verbose : bool
            Print progress
        
        Returns
        -------
        Dict
            Training results and metrics
            
        Raises
        ------
        ValidationError
            If parameters are invalid
        RuntimeError
            If training fails
        """
        # Validate inputs
        validate_positive_integer(n_datasets, 'n_datasets')
        if n_datasets < 10:
            raise ValidationError(
                f"n_datasets must be >= 10 for reliable training, got {n_datasets}"
            )
        
        valid_criteria = ['f1_score', 'auc', 'tpr', 'fdr_control']
        if criterion not in valid_criteria:
            raise ValidationError(
                f"criterion must be one of {valid_criteria}, got '{criterion}'"
            )
        
        from .simulation import simulate_diverse_datasets
        from .pipeline_runner import DEPipelineRunner
        from .profiler import RNAseqDataProfiler
        
        if verbose:
            print(f"\n{'='*60}")
            print("🦖 RAPTOR ML Training Pipeline")
            print(f"{'='*60}")
            print(f"Generating {n_datasets} diverse datasets...")
            print("This will take some time as we run real DE pipelines.\n")
        
        # Step 1: Generate diverse datasets
        logger.info(f"Generating {n_datasets} synthetic datasets...")
        try:
            datasets = simulate_diverse_datasets(n_datasets, seed=42)
        except Exception as e:
            raise RuntimeError(
                f"Failed to generate datasets: {e}\n"
                f"    Check simulation module is properly installed"
            )
        
        if not datasets:
            raise RuntimeError("No datasets generated")
        
        if verbose:
            print(f"✓ Generated {len(datasets)} datasets")
        
        # Step 2: Run pipelines and collect training data
        X_list = []
        y_list = []
        all_metrics = []
        
        runner = DEPipelineRunner()
        profiler = RNAseqDataProfiler()
        
        for i, dataset in enumerate(datasets):
            if verbose and (i + 1) % 10 == 0:
                print(f"  Processing dataset {i+1}/{len(datasets)}...")
            
            try:
                # Extract data
                counts = dataset['counts']
                metadata = dataset['metadata']
                true_de = dataset.get('true_de_genes', set())
                
                # Validate dataset
                validate_count_matrix(counts, allow_negative=False, min_genes=100, min_samples=4)
                validate_metadata(metadata, counts)
                
                # Profile the data
                profile = profiler.profile(counts, metadata)
                features = self._extract_features(profile)
                
                # Run all pipelines and evaluate
                results = runner.run_all_pipelines(counts, metadata)
                
                # Determine best pipeline based on criterion
                best_pipeline = self._determine_best_pipeline(
                    results, true_de, criterion
                )
                
                if best_pipeline:
                    X_list.append(features)
                    y_list.append(best_pipeline)
                    all_metrics.append(results)
            
            except Exception as e:
                logger.warning(f"Failed to process dataset {i+1}: {e}")
                continue
        
        if len(X_list) < 10:
            raise RuntimeError(
                f"Insufficient training data: only {len(X_list)} valid datasets. "
                f"Need at least 10 for reliable training."
            )
        
        if verbose:
            print(f"\n✓ Collected {len(X_list)} training samples")
        
        # Step 3: Prepare training data
        X = pd.DataFrame(X_list)
        y = np.array(y_list)
        
        # Check for missing features
        if X.isna().any().any():
            logger.warning("Missing values in features, filling with medians")
            X = X.fillna(X.median())
        
        self.feature_names = X.columns.tolist()
        
        # Encode labels
        y_encoded = self.label_encoder.fit_transform(y)
        
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y_encoded, test_size=0.2, stratify=y_encoded, random_state=42
        )
        
        # Scale features
        X_train_scaled = self.scaler.fit_transform(X_train)
        X_test_scaled = self.scaler.transform(X_test)
        
        # Train model
        if verbose:
            print(f"\nTraining {self.model_type} model...")
        
        self.model.fit(X_train_scaled, y_train)
        
        # Evaluate
        train_score = self.model.score(X_train_scaled, y_train)
        test_score = self.model.score(X_test_scaled, y_test)
        
        # Cross-validation
        cv_scores = cross_val_score(
            self.model, X_train_scaled, y_train,
            cv=min(5, len(np.unique(y_train))),
            n_jobs=-1
        )
        
        self.is_trained = True
        
        # Save training metadata
        self.training_metadata = {
            'n_datasets': n_datasets,
            'n_training_samples': len(X_list),
            'criterion': criterion,
            'model_type': self.model_type,
            'train_score': float(train_score),
            'test_score': float(test_score),
            'cv_mean': float(cv_scores.mean()),
            'cv_std': float(cv_scores.std()),
            'feature_names': self.feature_names,
            'pipelines': list(self.label_encoder.classes_),
            'trained_on': datetime.now().isoformat(),
            'version': self.VERSION
        }
        
        if verbose:
            print(f"\n{'='*60}")
            print("Training Complete!")
            print(f"{'='*60}")
            print(f"Train accuracy: {train_score:.3f}")
            print(f"Test accuracy:  {test_score:.3f}")
            print(f"CV score:       {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
            print(f"Pipelines:      {list(self.label_encoder.classes_)}")
        
        return self.training_metadata
    
    def _extract_features(self, profile: Union[Dict, 'DataProfile']) -> Dict[str, float]:
        """
        Extract ML features from data profile.
        
        Parameters
        ----------
        profile : Dict or DataProfile
            Data profile from RNAseqDataProfiler
            
        Returns
        -------
        Dict[str, float]
            Feature dictionary
            
        Raises
        ------
        ValidationError
            If profile is invalid
        """
        if not profile:
            raise ValidationError("Profile cannot be empty")
        
        # Convert DataProfile to dict if needed
        if hasattr(profile, '__dict__'):
            profile = profile.__dict__
        
        features = {}
        
        # Sample size features
        features['n_samples'] = profile.get('n_samples', 6)
        features['n_groups'] = profile.get('n_groups', 2)
        features['min_group_size'] = profile.get('min_group_size', 3)
        features['max_group_size'] = profile.get('max_group_size', 3)
        features['sample_size_balanced'] = float(
            profile.get('min_group_size', 3) >= profile.get('max_group_size', 3) * 0.8
        )
        
        # Dispersion features
        features['mean_bcv'] = profile.get('mean_bcv', 0.3)
        features['median_bcv'] = profile.get('median_bcv', 0.25)
        features['bcv_variance'] = profile.get('bcv_variance', 0.05)
        
        # Outlier features
        features['outlier_proportion'] = profile.get('outlier_proportion', 0.0)
        features['has_outliers'] = float(profile.get('has_outliers', False))
        
        # Library size features
        features['mean_library_size'] = profile.get('mean_library_size', 10_000_000)
        features['library_size_cv'] = profile.get('library_size_cv', 0.2)
        
        # Sparsity features
        features['sparsity'] = profile.get('sparsity', 0.5)
        features['zero_inflation_index'] = profile.get('zero_inflation_index', 0.0)
        
        # Count distribution features
        features['low_count_proportion'] = profile.get('low_count_proportion', 0.3)
        features['high_count_proportion'] = profile.get('high_count_proportion', 0.1)
        features['dynamic_range'] = profile.get('dynamic_range', 15.0)
        
        # Mean-variance relationship features
        features['mean_var_slope'] = profile.get('mean_var_slope', 1.5)
        features['poisson_fit_score'] = profile.get('poisson_fit_score', 0.5)
        
        # Validate all features are numeric
        for key, value in features.items():
            if not isinstance(value, (int, float)):
                raise ValidationError(
                    f"Feature '{key}' must be numeric, got {type(value)}"
                )
            if np.isnan(value) or np.isinf(value):
                raise ValidationError(
                    f"Feature '{key}' has invalid value: {value}"
                )
        
        return features
    
    def _determine_best_pipeline(
        self,
        results: Dict,
        true_de: set,
        criterion: str
    ) -> Optional[str]:
        """
        Determine best pipeline from results.
        
        Parameters
        ----------
        results : Dict
            Pipeline results
        true_de : set
            True DE gene set
        criterion : str
            Selection criterion
            
        Returns
        -------
        str or None
            Best pipeline name
        """
        pipeline_scores = {}
        
        for pipeline, metrics in results.items():
            if criterion == 'f1_score':
                score = metrics.get('f1_score', 0)
            elif criterion == 'auc':
                score = metrics.get('auc', 0.5)
            elif criterion == 'tpr':
                score = metrics.get('tpr', 0)
            elif criterion == 'fdr_control':
                score = 1 - metrics.get('fdr', 1)
            else:
                score = metrics.get('f1_score', 0)
            
            pipeline_scores[pipeline] = score
        
        if pipeline_scores:
            return max(pipeline_scores, key=pipeline_scores.get)
        return None
    
    @handle_errors(exit_on_error=False)
    def recommend(self, profile: Union[Dict, 'DataProfile']) -> MLRecommendation:
        """
        Get ML-based pipeline recommendation.
        
        Parameters
        ----------
        profile : Dict or DataProfile
            Data profile from RNAseqDataProfiler
        
        Returns
        -------
        MLRecommendation
            Recommendation with confidence and alternatives
            
        Raises
        ------
        RuntimeError
            If model not trained
        ValidationError
            If profile is invalid
        """
        if not self.is_trained:
            raise RuntimeError(
                "Model not trained. Call train() or load a trained model first."
            )
        
        if not profile:
            raise ValidationError("Profile cannot be empty")
        
        # Extract and validate features
        features = self._extract_features(profile)
        
        # Ensure correct feature order
        try:
            X = pd.DataFrame([features])[self.feature_names]
        except KeyError as e:
            raise ValidationError(
                f"Missing required features in profile: {e}\n"
                f"    Required features: {self.feature_names}"
            )
        
        # Check for missing values
        if X.isna().any().any():
            logger.warning("Missing values in features, filling with 0")
            X = X.fillna(0)
        
        # Scale and predict
        X_scaled = self.scaler.transform(X)
        proba = self.model.predict_proba(X_scaled)[0]
        classes = self.label_encoder.classes_
        
        # Sort by probability
        sorted_indices = np.argsort(proba)[::-1]
        
        best_idx = sorted_indices[0]
        best_pipeline = classes[best_idx]
        best_confidence = proba[best_idx]
        
        # Alternatives
        alternatives = [
            (classes[idx], proba[idx])
            for idx in sorted_indices[1:4]
        ]
        
        # Feature importances
        if hasattr(self.model, 'feature_importances_'):
            importances = dict(zip(
                self.feature_names,
                self.model.feature_importances_
            ))
            # Sort and take top
            importances = dict(sorted(
                importances.items(),
                key=lambda x: x[1],
                reverse=True
            )[:10])
        else:
            importances = {}
        
        return MLRecommendation(
            pipeline=best_pipeline,
            confidence=best_confidence,
            alternatives=alternatives,
            feature_importances=importances,
            model_version=self.VERSION
        )
    
    @handle_errors(exit_on_error=False)
    def save(self, path: str):
        """
        Save trained model.
        
        Parameters
        ----------
        path : str
            Output directory or file path
            
        Raises
        ------
        RuntimeError
            If model not trained
        PermissionError
            If cannot write to path
        """
        if not self.is_trained:
            raise RuntimeError("Cannot save untrained model. Train first with .train()")
        
        output_path = Path(path)
        
        # Validate output path is writable
        if output_path.suffix == '':
            # Directory
            validate_directory_path(output_path, must_exist=False, create_if_missing=True)
            model_file = output_path / f'ml_recommender_{self.model_type}.pkl'
        else:
            # File
            output_path.parent.mkdir(parents=True, exist_ok=True)
            model_file = output_path
        
        validate_output_writable(model_file)
        
        # Save everything in one file
        save_data = {
            'model': self.model,
            'scaler': self.scaler,
            'label_encoder': self.label_encoder,
            'feature_names': self.feature_names,
            'model_type': self.model_type,
            'metadata': self.training_metadata,
            'version': self.VERSION
        }
        
        try:
            joblib.dump(save_data, model_file)
            logger.info(f"✓ Model saved to {model_file}")
        except Exception as e:
            raise RuntimeError(f"Failed to save model: {e}")
        
        # Also save metadata as JSON
        meta_file = model_file.with_suffix('.json')
        try:
            with open(meta_file, 'w') as f:
                json.dump(self.training_metadata, f, indent=2)
            logger.info(f"✓ Metadata saved to {meta_file}")
        except Exception as e:
            logger.warning(f"Failed to save metadata JSON: {e}")
    
    @handle_errors(exit_on_error=False)
    def load(self, path: str):
        """
        Load trained model.
        
        Parameters
        ----------
        path : str
            Directory or file path
            
        Raises
        ------
        FileNotFoundError
            If model file not found
        RuntimeError
            If loading fails
        """
        input_path = Path(path)
        
        if input_path.is_dir():
            model_file = input_path / f'ml_recommender_{self.model_type}.pkl'
        else:
            model_file = input_path
        
        check_file_exists(model_file, "Model file")
        
        try:
            save_data = joblib.load(model_file)
        except Exception as e:
            raise RuntimeError(
                f"Failed to load model from {model_file}: {e}\n"
                f"    File may be corrupted or incompatible"
            )
        
        # Validate loaded data
        required_keys = ['model', 'scaler', 'label_encoder', 'feature_names']
        missing_keys = set(required_keys) - set(save_data.keys())
        if missing_keys:
            raise RuntimeError(
                f"Invalid model file: missing keys {missing_keys}\n"
                f"    File may be from an older version"
            )
        
        self.model = save_data['model']
        self.scaler = save_data['scaler']
        self.label_encoder = save_data['label_encoder']
        self.feature_names = save_data['feature_names']
        self.model_type = save_data.get('model_type', 'random_forest')
        self.training_metadata = save_data.get('metadata', {})
        self.is_trained = True
        
        logger.info(f"✓ Model loaded from {model_file}")
        logger.info(f"  Model type: {self.model_type}")
        logger.info(f"  Trained on: {self.training_metadata.get('trained_on', 'unknown')}")
        logger.info(f"  Pipelines: {self.training_metadata.get('pipelines', [])}")


# =============================================================================
# Convenience Functions
# =============================================================================

@handle_errors(exit_on_error=False)
def train_ml_recommender(
    n_datasets: int = 100,
    model_type: str = 'random_forest',
    output_path: str = 'models/ml_recommender.pkl',
    verbose: bool = True
) -> MLPipelineRecommender:
    """
    Convenience function to train and save ML recommender.
    
    Parameters
    ----------
    n_datasets : int
        Number of training datasets (must be >= 10)
    model_type : str
        'random_forest' or 'gradient_boosting'
    output_path : str
        Where to save model
    verbose : bool
        Print progress
    
    Returns
    -------
    MLPipelineRecommender
        Trained recommender
        
    Raises
    ------
    ValidationError
        If parameters invalid
    RuntimeError
        If training fails
    """
    # Validate inputs
    validate_positive_integer(n_datasets, 'n_datasets')
    if n_datasets < 10:
        raise ValidationError(
            f"n_datasets must be >= 10, got {n_datasets}"
        )
    
    valid_models = ['random_forest', 'gradient_boosting']
    if model_type not in valid_models:
        raise ValidationError(
            f"model_type must be one of {valid_models}, got '{model_type}'"
        )
    
    # Validate output path
    validate_output_writable(Path(output_path))
    
    # Train
    recommender = MLPipelineRecommender(model_type=model_type)
    results = recommender.train(n_datasets=n_datasets, verbose=verbose)
    recommender.save(output_path)
    
    return recommender


@handle_errors(exit_on_error=False)
def quick_ml_recommend(
    profile: Union[Dict, 'DataProfile'],
    model_path: str
) -> MLRecommendation:
    """
    Quick ML recommendation from saved model.
    
    Parameters
    ----------
    profile : Dict or DataProfile
        Data profile
    model_path : str
        Path to trained model
    
    Returns
    -------
    MLRecommendation
        Pipeline recommendation
        
    Raises
    ------
    FileNotFoundError
        If model not found
    ValidationError
        If profile invalid
    """
    if not profile:
        raise ValidationError("Profile cannot be empty")
    
    # Validate model path
    path = Path(model_path)
    if not path.exists():
        raise FileNotFoundError(
            f"Model not found: {model_path}\n"
            f"    Train a model first with train_ml_recommender()"
        )
    
    recommender = MLPipelineRecommender(model_path=model_path)
    return recommender.recommend(profile)


if __name__ == '__main__':
    print("🦖 RAPTOR ML-Based Pipeline Recommender (VALIDATED)")
    print("=" * 60)
    print("\nTRUE Machine Learning - trains on actual benchmark results!")
    print("With comprehensive input validation and error handling.")
    print("\nUsage:")
    print("  # Train new model (takes time - runs real pipelines)")
    print("  from ml_recommender import MLPipelineRecommender")
    print("  recommender = MLPipelineRecommender()")
    print("  recommender.train(n_datasets=100)")
    print("  recommender.save('models/ml_recommender.pkl')")
    print("")
    print("  # Use trained model")
    print("  recommender = MLPipelineRecommender(model_path='models/')")
    print("  rec = recommender.recommend(profile)")
    print("  print(rec.pipeline, rec.confidence)")
