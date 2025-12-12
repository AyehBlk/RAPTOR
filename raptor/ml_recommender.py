"""
ML-Based Pipeline Recommender

Machine learning-powered pipeline recommendation system that learns from
historical benchmark data to predict optimal pipelines for new datasets.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import logging
import json
import pickle
from datetime import datetime

# ML libraries
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score, train_test_split
from sklearn.metrics import classification_report, confusion_matrix
import joblib

logger = logging.getLogger(__name__)


class FeatureExtractor:
    """
    Extract and engineer features from RNA-seq data profiles.
    
    Converts raw data characteristics into ML-ready features for
    pipeline recommendation.
    """
    
    @staticmethod
    def extract_features(profile: Dict) -> pd.DataFrame:
        """
        Extract feature vector from data profile.
        
        Parameters
        ----------
        profile : dict
            Data profile from RNAseqDataProfiler
        
        Returns
        -------
        pd.DataFrame
            Feature vector (single row)
        """
        features = {}
        
        # Design features
        if 'design' in profile:
            features['n_samples'] = profile['design'].get('n_samples', 0)
            features['n_genes'] = profile['design'].get('n_genes', 0)
            features['n_conditions'] = profile['design'].get('n_conditions', 2)
            features['samples_per_condition'] = profile['design'].get('samples_per_condition', 3)
            features['has_replicates'] = int(features['samples_per_condition'] > 1)
            features['is_paired'] = int(profile['design'].get('is_paired', False))
            
        # Library statistics
        if 'library_stats' in profile:
            features['mean_lib_size'] = profile['library_stats'].get('mean', 0)
            features['median_lib_size'] = profile['library_stats'].get('median', 0)
            features['lib_size_cv'] = profile['library_stats'].get('cv', 0)
            features['lib_size_range'] = profile['library_stats'].get('range', 0)
            features['lib_size_skewness'] = profile['library_stats'].get('skewness', 0)
        
        # Count distribution
        if 'count_distribution' in profile:
            features['zero_pct'] = profile['count_distribution'].get('zero_pct', 0)
            features['low_count_pct'] = profile['count_distribution'].get('low_count_pct', 0)
            features['mean_count'] = profile['count_distribution'].get('mean', 0)
            features['median_count'] = profile['count_distribution'].get('median', 0)
            features['count_variance'] = profile['count_distribution'].get('variance', 0)
        
        # Expression distribution
        if 'expression_distribution' in profile:
            features['high_expr_genes'] = profile['expression_distribution'].get('high_expr_genes', 0)
            features['medium_expr_genes'] = profile['expression_distribution'].get('medium_expr_genes', 0)
            features['low_expr_genes'] = profile['expression_distribution'].get('low_expr_genes', 0)
            features['expr_dynamic_range'] = profile['expression_distribution'].get('dynamic_range', 0)
        
        # Biological variation
        if 'biological_variation' in profile:
            features['bcv'] = profile['biological_variation'].get('bcv', 0)
            features['dispersion_mean'] = profile['biological_variation'].get('dispersion_mean', 0)
            features['dispersion_trend'] = profile['biological_variation'].get('dispersion_trend', 0)
            features['outlier_genes'] = profile['biological_variation'].get('outlier_genes', 0)
        
        # Sequencing depth
        if 'sequencing' in profile:
            features['total_reads'] = profile['sequencing'].get('total_reads', 0)
            features['reads_per_gene'] = profile['sequencing'].get('reads_per_gene', 0)
            depth_map = {'low': 1, 'medium': 2, 'high': 3, 'very_high': 4}
            features['depth_category'] = depth_map.get(
                profile['sequencing'].get('depth_category', 'medium'), 2
            )
        
        # Data complexity features
        if 'complexity' in profile:
            features['complexity_score'] = profile['complexity'].get('score', 50)
            features['noise_level'] = profile['complexity'].get('noise_level', 0)
            features['signal_strength'] = profile['complexity'].get('signal_strength', 0)
        
        # Derived features
        features['samples_genes_ratio'] = features.get('n_samples', 1) / max(features.get('n_genes', 1), 1)
        features['effective_lib_size'] = features.get('mean_lib_size', 0) * (1 - features.get('zero_pct', 0) / 100)
        features['data_quality_score'] = FeatureExtractor._calculate_quality_score(features)
        
        return pd.DataFrame([features])
    
    @staticmethod
    def _calculate_quality_score(features: Dict) -> float:
        """
        Calculate overall data quality score (0-100).
        
        Parameters
        ----------
        features : dict
            Extracted features
        
        Returns
        -------
        float
            Quality score
        """
        score = 100.0
        
        # Penalize high library size CV
        lib_cv = features.get('lib_size_cv', 0)
        if lib_cv > 0.3:
            score -= min(20, (lib_cv - 0.3) * 50)
        
        # Penalize high zero percentage
        zero_pct = features.get('zero_pct', 0)
        if zero_pct > 60:
            score -= min(25, (zero_pct - 60) / 2)
        
        # Penalize high BCV
        bcv = features.get('bcv', 0)
        if bcv > 0.4:
            score -= min(15, (bcv - 0.4) * 30)
        
        # Bonus for good sample size
        n_samples = features.get('n_samples', 0)
        if n_samples >= 6:
            score += min(10, n_samples)
        
        return max(0, min(100, score))


class MLPipelineRecommender:
    """
    Machine learning-based pipeline recommender.
    
    Uses historical benchmark data to train models that predict the best
    pipeline for new datasets based on their characteristics.
    
    Parameters
    ----------
    model_type : str
        'random_forest' or 'gradient_boosting'
    model_path : str, optional
        Path to pre-trained model
    
    Examples
    --------
    >>> recommender = MLPipelineRecommender(model_type='random_forest')
    >>> recommender.train_from_benchmarks('benchmark_history/')
    >>> recommendation = recommender.recommend(data_profile)
    """
    
    def __init__(
        self,
        model_type: str = 'random_forest',
        model_path: Optional[str] = None
    ):
        """Initialize ML recommender."""
        self.model_type = model_type
        self.model = None
        self.scaler = StandardScaler()
        self.feature_names = None
        self.training_history = []
        
        if model_path and Path(model_path).exists():
            self.load_model(model_path)
        else:
            self._initialize_model()
        
        logger.info(f"Initialized ML recommender: {model_type}")
    
    def _initialize_model(self):
        """Initialize ML model."""
        if self.model_type == 'random_forest':
            self.model = RandomForestClassifier(
                n_estimators=200,
                max_depth=15,
                min_samples_split=5,
                min_samples_leaf=2,
                max_features='sqrt',
                random_state=42,
                n_jobs=-1,
                class_weight='balanced'
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
            raise ValueError(f"Unknown model type: {self.model_type}")
    
    def train_from_benchmarks(
        self,
        benchmark_dir: str,
        performance_metric: str = 'f1_score'
    ) -> Dict:
        """
        Train model from historical benchmark results.
        
        Parameters
        ----------
        benchmark_dir : str
            Directory containing benchmark results
        performance_metric : str
            Metric to optimize: 'f1_score', 'accuracy', 'runtime', or 'combined'
        
        Returns
        -------
        dict
            Training results and performance metrics
        """
        logger.info(f"Training from benchmarks: {benchmark_dir}")
        
        # Load all benchmark results
        X, y, metadata = self._load_benchmark_data(benchmark_dir, performance_metric)
        
        if len(X) < 10:
            raise ValueError(
                f"Insufficient training data: {len(X)} samples. "
                "Need at least 10 benchmark results."
            )
        
        logger.info(f"Loaded {len(X)} benchmark results")
        logger.info(f"Pipeline distribution: {pd.Series(y).value_counts().to_dict()}")
        
        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=0.2, random_state=42, stratify=y
        )
        
        # Scale features
        X_train_scaled = self.scaler.fit_transform(X_train)
        X_test_scaled = self.scaler.transform(X_test)
        
        # Train model
        logger.info("Training model...")
        self.model.fit(X_train_scaled, y_train)
        
        # Evaluate
        train_score = self.model.score(X_train_scaled, y_train)
        test_score = self.model.score(X_test_scaled, y_test)
        
        # Cross-validation
        cv_scores = cross_val_score(
            self.model, X_train_scaled, y_train, cv=5, n_jobs=-1
        )
        
        # Predictions
        y_pred = self.model.predict(X_test_scaled)
        
        # Feature importance
        if hasattr(self.model, 'feature_importances_'):
            feature_importance = pd.DataFrame({
                'feature': self.feature_names,
                'importance': self.model.feature_importances_
            }).sort_values('importance', ascending=False)
        else:
            feature_importance = None
        
        results = {
            'model_type': self.model_type,
            'n_samples': len(X),
            'n_features': X.shape[1],
            'train_score': float(train_score),
            'test_score': float(test_score),
            'cv_mean': float(cv_scores.mean()),
            'cv_std': float(cv_scores.std()),
            'classification_report': classification_report(y_test, y_pred, output_dict=True),
            'confusion_matrix': confusion_matrix(y_test, y_pred).tolist(),
            'feature_importance': feature_importance.to_dict('records') if feature_importance is not None else None,
            'timestamp': datetime.now().isoformat()
        }
        
        self.training_history.append(results)
        
        logger.info(f"Training complete - Test accuracy: {test_score:.3f}")
        logger.info(f"Cross-validation: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
        
        return results
    
    def _load_benchmark_data(
        self,
        benchmark_dir: str,
        performance_metric: str
    ) -> Tuple[pd.DataFrame, np.ndarray, List[Dict]]:
        """
        Load and process benchmark data.
        
        Parameters
        ----------
        benchmark_dir : str
            Directory with benchmark results
        performance_metric : str
            Metric to use for determining best pipeline
        
        Returns
        -------
        tuple
            (features_df, labels_array, metadata_list)
        """
        benchmark_path = Path(benchmark_dir)
        
        X_list = []
        y_list = []
        metadata_list = []
        
        # Find all benchmark result files
        for result_file in benchmark_path.glob('**/benchmark_results.json'):
            try:
                with open(result_file, 'r') as f:
                    benchmark_data = json.load(f)
                
                # Load corresponding profile
                profile_file = result_file.parent / 'data_profile.json'
                if not profile_file.exists():
                    logger.warning(f"Profile not found for {result_file}")
                    continue
                
                with open(profile_file, 'r') as f:
                    profile = json.load(f)
                
                # Extract features
                features = FeatureExtractor.extract_features(profile)
                
                # Determine best pipeline
                best_pipeline = self._determine_best_pipeline(
                    benchmark_data,
                    performance_metric
                )
                
                if best_pipeline is not None:
                    X_list.append(features)
                    y_list.append(best_pipeline)
                    metadata_list.append({
                        'source': str(result_file),
                        'profile': profile.get('summary', {}),
                        'best_pipeline': best_pipeline
                    })
                
            except Exception as e:
                logger.warning(f"Error loading {result_file}: {e}")
                continue
        
        if not X_list:
            raise ValueError(f"No valid benchmark data found in {benchmark_dir}")
        
        # Combine features
        X = pd.concat(X_list, ignore_index=True)
        y = np.array(y_list)
        
        # Store feature names
        self.feature_names = X.columns.tolist()
        
        return X, y, metadata_list
    
    def _determine_best_pipeline(
        self,
        benchmark_data: Dict,
        metric: str
    ) -> Optional[int]:
        """
        Determine best pipeline from benchmark results.
        
        Parameters
        ----------
        benchmark_data : dict
            Benchmark results
        metric : str
            Performance metric to use
        
        Returns
        -------
        int or None
            Best pipeline ID
        """
        scores = {}
        
        for pipeline_id, result in benchmark_data.items():
            if isinstance(pipeline_id, str):
                pipeline_id = int(pipeline_id)
            
            if result.get('status') != 'success':
                continue
            
            # Calculate score based on metric
            if metric == 'f1_score':
                score = result.get('metrics', {}).get('f1_score', 0)
            elif metric == 'accuracy':
                score = result.get('metrics', {}).get('accuracy', 0)
            elif metric == 'runtime':
                # Lower is better for runtime
                runtime = result.get('runtime', float('inf'))
                score = 1.0 / (runtime + 1)  # Avoid division by zero
            elif metric == 'combined':
                # Weighted combination
                f1 = result.get('metrics', {}).get('f1_score', 0)
                runtime = result.get('runtime', float('inf'))
                score = f1 * 0.7 + (1.0 / (runtime + 1)) * 0.3
            else:
                score = result.get('score', 0)
            
            scores[pipeline_id] = score
        
        if not scores:
            return None
        
        return max(scores.items(), key=lambda x: x[1])[0]
    
    def recommend(self, profile: Dict, top_k: int = 3) -> Dict:
        """
        Recommend pipeline for given data profile.
        
        Parameters
        ----------
        profile : dict
            Data profile from RNAseqDataProfiler
        top_k : int
            Number of top recommendations to return
        
        Returns
        -------
        dict
            Recommendation with confidence scores
        """
        if self.model is None:
            raise ValueError("Model not trained. Call train_from_benchmarks() first.")
        
        # Extract features
        features = FeatureExtractor.extract_features(profile)
        
        # Ensure all expected features are present
        for col in self.feature_names:
            if col not in features.columns:
                features[col] = 0
        
        # Reorder to match training
        features = features[self.feature_names]
        
        # Scale
        features_scaled = self.scaler.transform(features)
        
        # Predict probabilities
        probabilities = self.model.predict_proba(features_scaled)[0]
        classes = self.model.classes_
        
        # Get top k recommendations
        top_indices = np.argsort(probabilities)[-top_k:][::-1]
        
        recommendations = []
        for idx in top_indices:
            pipeline_id = int(classes[idx])
            confidence = float(probabilities[idx])
            
            recommendations.append({
                'pipeline_id': pipeline_id,
                'pipeline_name': self._get_pipeline_name(pipeline_id),
                'confidence': confidence,
                'reasons': self._generate_reasons(pipeline_id, profile, confidence)
            })
        
        # Primary recommendation
        primary = recommendations[0]
        
        result = {
            'pipeline_id': primary['pipeline_id'],
            'pipeline_name': primary['pipeline_name'],
            'confidence': primary['confidence'],
            'model_type': self.model_type,
            'reasons': primary['reasons'],
            'alternatives': recommendations[1:] if len(recommendations) > 1 else [],
            'feature_contributions': self._get_feature_contributions(features_scaled[0]),
            'timestamp': datetime.now().isoformat()
        }
        
        logger.info(
            f"Recommended Pipeline {primary['pipeline_id']} "
            f"with {primary['confidence']:.1%} confidence"
        )
        
        return result
    
    def _get_feature_contributions(self, features_scaled: np.ndarray) -> List[Dict]:
        """
        Get feature contributions to prediction.
        
        Parameters
        ----------
        features_scaled : np.ndarray
            Scaled feature vector
        
        Returns
        -------
        list of dict
            Top contributing features
        """
        if not hasattr(self.model, 'feature_importances_'):
            return []
        
        # Get importances
        importances = self.model.feature_importances_
        
        # Weight by feature values
        contributions = importances * np.abs(features_scaled)
        
        # Get top features
        top_indices = np.argsort(contributions)[-5:][::-1]
        
        result = []
        for idx in top_indices:
            result.append({
                'feature': self.feature_names[idx],
                'importance': float(importances[idx]),
                'value': float(features_scaled[idx]),
                'contribution': float(contributions[idx])
            })
        
        return result
    
    def _generate_reasons(
        self,
        pipeline_id: int,
        profile: Dict,
        confidence: float
    ) -> List[str]:
        """Generate human-readable reasons for recommendation."""
        reasons = []
        
        # Add confidence statement
        if confidence > 0.7:
            reasons.append(f"High confidence recommendation ({confidence:.1%})")
        elif confidence > 0.5:
            reasons.append(f"Moderate confidence recommendation ({confidence:.1%})")
        else:
            reasons.append(f"Low confidence - consider alternatives ({confidence:.1%})")
        
        # Pipeline-specific reasons based on profile
        n_samples = profile.get('design', {}).get('n_samples', 0)
        bcv = profile.get('biological_variation', {}).get('bcv', 0)
        depth = profile.get('sequencing', {}).get('depth_category', 'medium')
        
        if pipeline_id == 1:  # STAR-RSEM-DESeq2
            reasons.append("Gold standard pipeline with excellent accuracy")
            if n_samples >= 6:
                reasons.append("Well-suited for your sample size")
        
        elif pipeline_id == 3:  # Salmon-edgeR
            reasons.append("Fast pseudo-alignment with robust statistics")
            if depth in ['high', 'very_high']:
                reasons.append("Efficient for high-depth sequencing")
        
        elif pipeline_id == 4:  # Kallisto-Sleuth
            reasons.append("Ultra-fast quantification with uncertainty modeling")
            if n_samples >= 6:
                reasons.append("Sleuth benefits from multiple replicates")
        
        elif pipeline_id == 6:  # NOISeq
            reasons.append("Designed for datasets with limited replicates")
            if n_samples < 3:
                reasons.append("Optimal for your low replicate count")
        
        # Add data-specific reasons
        if bcv > 0.4:
            reasons.append("Handles high biological variability well")
        
        return reasons
    
    @staticmethod
    def _get_pipeline_name(pipeline_id: int) -> str:
        """Get pipeline name from ID."""
        names = {
            1: "STAR-RSEM-DESeq2",
            2: "HISAT2-StringTie-Ballgown",
            3: "Salmon-edgeR",
            4: "Kallisto-Sleuth",
            5: "STAR-HTSeq-limma-voom",
            6: "STAR-featureCounts-NOISeq",
            7: "Bowtie2-RSEM-EBSeq",
            8: "HISAT2-Cufflinks-Cuffdiff"
        }
        return names.get(pipeline_id, f"Pipeline {pipeline_id}")
    
    def save_model(self, output_dir: str):
        """
        Save trained model and scaler.
        
        Parameters
        ----------
        output_dir : str
            Output directory
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Save model
        model_file = output_path / f'ml_recommender_{self.model_type}.pkl'
        joblib.dump(self.model, model_file)
        
        # Save scaler
        scaler_file = output_path / f'scaler_{self.model_type}.pkl'
        joblib.dump(self.scaler, scaler_file)
        
        # Save metadata
        metadata = {
            'model_type': self.model_type,
            'feature_names': self.feature_names,
            'training_history': self.training_history,
            'timestamp': datetime.now().isoformat()
        }
        
        metadata_file = output_path / f'metadata_{self.model_type}.json'
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        logger.info(f"Model saved to {output_path}")
    
    def load_model(self, model_dir: str):
        """
        Load pre-trained model.
        
        Parameters
        ----------
        model_dir : str
            Directory containing model files
        """
        model_path = Path(model_dir)
        
        # Load model
        model_file = model_path / f'ml_recommender_{self.model_type}.pkl'
        if not model_file.exists():
            raise FileNotFoundError(f"Model not found: {model_file}")
        
        self.model = joblib.load(model_file)
        
        # Load scaler
        scaler_file = model_path / f'scaler_{self.model_type}.pkl'
        self.scaler = joblib.load(scaler_file)
        
        # Load metadata
        metadata_file = model_path / f'metadata_{self.model_type}.json'
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)
        
        self.feature_names = metadata['feature_names']
        self.training_history = metadata.get('training_history', [])
        
        logger.info(f"Model loaded from {model_path}")
    
    def evaluate_on_new_data(self, benchmark_dir: str, metric: str = 'f1_score') -> Dict:
        """
        Evaluate model on new benchmark data.
        
        Parameters
        ----------
        benchmark_dir : str
            Directory with new benchmark results
        metric : str
            Performance metric
        
        Returns
        -------
        dict
            Evaluation results
        """
        X, y_true, metadata = self._load_benchmark_data(benchmark_dir, metric)
        
        # Predict
        X_scaled = self.scaler.transform(X)
        y_pred = self.model.predict(X_scaled)
        
        # Calculate metrics
        from sklearn.metrics import accuracy_score, precision_recall_fscore_support
        
        accuracy = accuracy_score(y_true, y_pred)
        precision, recall, f1, _ = precision_recall_fscore_support(
            y_true, y_pred, average='weighted'
        )
        
        results = {
            'n_samples': len(X),
            'accuracy': float(accuracy),
            'precision': float(precision),
            'recall': float(recall),
            'f1_score': float(f1),
            'confusion_matrix': confusion_matrix(y_true, y_pred).tolist(),
            'timestamp': datetime.now().isoformat()
        }
        
        logger.info(f"Evaluation on new data - Accuracy: {accuracy:.3f}")
        
        return results


# =============================================================================
# Convenience Functions
# =============================================================================

def train_recommender(
    benchmark_dir: str,
    model_type: str = 'random_forest',
    output_dir: str = 'models/'
) -> MLPipelineRecommender:
    """
    Convenience function to train and save recommender.
    
    Parameters
    ----------
    benchmark_dir : str
        Directory with benchmark results
    model_type : str
        'random_forest' or 'gradient_boosting'
    output_dir : str
        Where to save model
    
    Returns
    -------
    MLPipelineRecommender
        Trained recommender
    """
    recommender = MLPipelineRecommender(model_type=model_type)
    results = recommender.train_from_benchmarks(benchmark_dir)
    
    print("\n=== Training Results ===")
    print(f"Model: {model_type}")
    print(f"Samples: {results['n_samples']}")
    print(f"Test Accuracy: {results['test_score']:.3f}")
    print(f"CV Score: {results['cv_mean']:.3f} ± {results['cv_std']:.3f}")
    
    recommender.save_model(output_dir)
    
    return recommender


if __name__ == '__main__':
    print("RAPTOR ML-Based Pipeline Recommender")
    print("====================================")
    print("\nMachine learning-powered pipeline recommendation.")
    print("\nFeatures:")
    print("  • Learns from historical benchmark data")
    print("  • Predicts optimal pipeline with confidence scores")
    print("  • Provides interpretable recommendations")
    print("  • Supports Random Forest and Gradient Boosting")
    print("\nUsage:")
    print("  from raptor.ml_recommender import MLPipelineRecommender")
    print("  recommender = MLPipelineRecommender()")
    print("  recommender.train_from_benchmarks('benchmark_history/')")
    print("  recommendation = recommender.recommend(data_profile)")
