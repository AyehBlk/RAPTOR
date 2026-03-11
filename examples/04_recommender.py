#!/usr/bin/env python3
"""
RAPTOR v2.2.0 Example Script: Pipeline Recommender (Module 4)

Demonstrates pipeline recommendation based on data characteristics.
Supports both rule-based and ML-based recommendations.

This is Module 4 of the RAPTOR workflow (Stage 1: Fast Profiling):
  M1: Quantify (FASTQ → quick_gene_counts.csv)
  M2: Sample QC (Quality Assessment & Outlier Detection)
  M3: Profile (Data Profiling - 32 features)
  M4: Recommend (Pipeline Recommendation) ← THIS SCRIPT

Recommends one of these DE analysis pipelines:
  - DESeq2: General purpose, handles batch effects
  - edgeR: Handles low counts and high dispersion
  - limma-voom: Best for large samples and complex designs
  - Wilcoxon: Non-parametric for large samples
  - edgeR_robust: Handles outliers well

Input: results/profile.json (or counts file for auto-profiling)
Output: results/recommendation.json

After recommendation, run the suggested pipeline in R for DE analysis.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
License: MIT
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

# =============================================================================
# CONSTANTS - Architecture Compliant (v2.2.0)
# =============================================================================
DEFAULT_CLEAN_COUNTS = "results/qc/counts_clean.csv"
DEFAULT_PROFILE = "results/profile.json"
DEFAULT_OUTPUT_DIR = "results"

# Check for dependencies
try:
    import numpy as np
    import pandas as pd
except ImportError:
    print("ERROR: numpy and pandas are required")
    print("Install with: pip install numpy pandas")
    sys.exit(1)

# RAPTOR imports - UPDATED FOR v2.2.0
RAPTOR_AVAILABLE = True
try:
    from raptor import (
        # Recommender classes
        PipelineRecommender,
        Recommendation,
        
        # Convenience function
        recommend_pipeline,
        
        # For auto-profiling
        RNAseqDataProfiler,
        DataProfile,
        
        # Validation
        validate_count_matrix,
        validate_file_path,
        validate_directory_path,
        ValidationError
    )
except ImportError:
    RAPTOR_AVAILABLE = False
    print("NOTE: RAPTOR not installed. Running in demo mode only.")
    print("Install RAPTOR with: pip install -e .")
    
    # Create dummy classes for demo
    class PipelineRecommender:
        def __init__(self, profile): 
            self.profile = profile
        def get_recommendation(self): 
            return None
    
    class Recommendation:
        def __init__(self, **kwargs):
            for k, v in kwargs.items():
                setattr(self, k, v)
        def summary(self):
            return "Demo recommendation summary"
    
    class RNAseqDataProfiler:
        def __init__(self, counts, metadata=None, group_column=None): 
            pass
        def run_full_profile(self): 
            return None
    
    class DataProfile:
        pass
    
    def recommend_pipeline(profile): 
        return None
    def validate_count_matrix(df, **kwargs): 
        pass
    def validate_file_path(p, **kwargs): 
        return Path(p)
    def validate_directory_path(p, **kwargs): 
        return Path(p)
    class ValidationError(ValueError): 
        pass


def print_banner():
    """Print RAPTOR banner."""
    print("""
    ╔══════════════════════════════════════════════════════════════╗
    ║       🦖 RAPTOR v2.2.0 - Pipeline Recommender (Module 4)     ║
    ║                                                              ║
    ║   Rule-Based Pipeline Recommendation for DE Analysis        ║
    ║   Recommends: DESeq2 | edgeR | limma-voom | Wilcoxon        ║
    ╚══════════════════════════════════════════════════════════════╝
    """)


def display_recommendation(recommendation):
    """Display recommendation summary."""
    print("\n🎯 Pipeline Recommendation:")
    print("  " + "="*70)
    
    if recommendation is None:
        # Demo mode
        print("  🥇 PRIMARY RECOMMENDATION: DESeq2")
        print("     Confidence: 88%")
        print("     Reason: DESeq2: appropriate for small samples, suitable for moderate dispersion")
        
        print("\n  🥈 ALTERNATIVE: edgeR")
        print("     Confidence: 82%")
        print("     Reason: edgeR: fast and flexible NB model")
        
        print("\n  Decision Factors:")
        print("     • Sample size per group: 6")
        print("     • BCV (biological variation): 0.32 (moderate)")
        print("     • Sparsity: 0.15 (15% zeros)")
        
        print("  " + "="*70)
        return
    
    # Real recommendation - use the summary() method if available
    if hasattr(recommendation, 'summary'):
        summary_text = recommendation.summary()
        print(summary_text)
    else:
        # Fallback display
        primary = getattr(recommendation, 'primary_pipeline', 'DESeq2')
        score = getattr(recommendation, 'primary_score', 88.0)
        reason = getattr(recommendation, 'primary_reason', '')
        
        print(f"  🥇 PRIMARY RECOMMENDATION: {primary}")
        print(f"     Confidence: {score:.0f}%")
        if reason:
            print(f"     Reason: {reason}")
        
        alt = getattr(recommendation, 'alternative_pipeline', None)
        if alt:
            alt_score = getattr(recommendation, 'alternative_score', 0)
            alt_reason = getattr(recommendation, 'alternative_reason', '')
            print(f"\n  🥈 ALTERNATIVE: {alt}")
            print(f"     Confidence: {alt_score:.0f}%")
            if alt_reason:
                print(f"     Reason: {alt_reason}")
        
        # Show warnings if any
        warnings = getattr(recommendation, 'warnings', [])
        if warnings:
            print(f"\n  ⚠️  WARNINGS:")
            for warning in warnings[:3]:
                print(f"     • {warning}")
        
        print("  " + "="*70)


def run_recommendation(profile_file=None, counts_file=None, metadata_file=None,
                      output_dir=None, demo=False, verbose=False):
    """
    Run pipeline recommendation with validation.
    
    Parameters
    ----------
    profile_file : str, optional
        Path to pre-computed profile JSON
    counts_file : str, optional
        Path to count matrix for auto-profiling
    metadata_file : str, optional
        Path to metadata for auto-profiling
    output_dir : str, optional
        Output directory
    demo : bool
        Run in demo mode
    verbose : bool
        Show detailed reasoning
    """
    
    profile = None
    
    # Mode 1: Load existing profile
    if profile_file and not demo:
        try:
            print("\n🔍 Loading profile...")
            
            profile_path = validate_file_path(profile_file, must_exist=True)
            print(f"  ✓ Profile file found: {profile_path}")
            
            with open(profile_path, 'r') as f:
                profile_data = json.load(f)
            
            print(f"  ✓ Profile loaded")
            
            # Convert to DataProfile if RAPTOR available
            if RAPTOR_AVAILABLE:
                # Check if it's already in the right format
                if 'features' in profile_data:
                    profile_dict = profile_data.get('features', profile_data)
                else:
                    profile_dict = profile_data
                
                # Create DataProfile object with required fields
                profile = type('Profile', (), profile_dict)()
            else:
                profile = type('Profile', (), profile_data)()
            
        except FileNotFoundError:
            print(f"❌ Profile file not found: {profile_file}")
            print(f"   Run Module 3 first: python 03_data_profiler.py")
            sys.exit(1)
        except json.JSONDecodeError as e:
            print(f"❌ Invalid JSON in profile file: {e}")
            sys.exit(1)
        except Exception as e:
            print(f"❌ Error loading profile: {e}")
            sys.exit(1)
    
    # Mode 2: Profile from counts
    elif counts_file and not demo:
        try:
            print("\n🔍 Loading counts for auto-profiling...")
            
            counts_path = validate_file_path(counts_file, must_exist=True)
            print(f"  ✓ Count file found: {counts_path}")
            
            counts_df = pd.read_csv(counts_path, index_col=0)
            print(f"  ✓ Count matrix loaded: {counts_df.shape[0]} genes × {counts_df.shape[1]} samples")
            
            validate_count_matrix(counts_df, min_genes=10, min_samples=2)
            print(f"  ✓ Count matrix validated")
            
            if not RAPTOR_AVAILABLE:
                print("❌ RAPTOR required for auto-profiling")
                print("   Install with: pip install -e .")
                sys.exit(1)
            
            # Load metadata if provided
            metadata_df = None
            if metadata_file:
                try:
                    metadata_df = pd.read_csv(metadata_file)
                    print(f"  ✓ Metadata loaded: {len(metadata_df)} samples")
                except Exception as e:
                    print(f"  ⚠️  Warning: Could not load metadata: {e}")
            
            print("\n📊 Running auto-profiling...")
            profiler = RNAseqDataProfiler(counts_df, metadata_df)
            profile = profiler.run_full_profile()
            print("  ✓ Profile generated")
            
        except FileNotFoundError:
            print(f"❌ Count file not found: {counts_file}")
            sys.exit(1)
        except ValidationError as e:
            print(f"❌ Count matrix validation failed: {e}")
            sys.exit(1)
        except Exception as e:
            print(f"❌ Error during auto-profiling: {e}")
            sys.exit(1)
    
    # Mode 3: Demo mode
    else:
        print("\n🎬 Running in DEMO mode")
        # Create demo profile with realistic values
        profile = type('Profile', (), {
            'n_samples': 12,
            'n_genes': 20000,
            'n_groups': 2,
            'min_group_size': 6,
            'bcv': 0.32,
            'bcv_category': 'moderate',
            'library_size_cv': 0.12,
            'sparsity': 0.15,
            'low_count_proportion': 0.25,
            'has_outliers': False,
            'has_batch_effect': False,
            'design_complexity': 'simple'
        })()
    
    # Validate output directory
    try:
        output_path = validate_directory_path(
            output_dir or DEFAULT_OUTPUT_DIR,
            create_if_missing=True
        )
        print(f"  ✓ Output directory: {output_path}")
    except Exception as e:
        print(f"❌ Error with output directory: {e}")
        sys.exit(1)
    
    # Run recommendation
    print(f"\n🎯 Generating recommendation...")
    
    if not RAPTOR_AVAILABLE or demo:
        print("  (Demo mode - showing simulated recommendation)")
        
        # Demo recommendation
        recommendation = Recommendation(
            primary_pipeline='DESeq2',
            primary_score=88.0,
            primary_reason='DESeq2: appropriate for small samples, suitable for moderate dispersion',
            alternative_pipeline='edgeR',
            alternative_score=82.0,
            alternative_reason='edgeR: fast and flexible NB model',
            third_option='limma-voom',
            third_score=75.0,
            decision_factors={
                'sample_size_per_group': 6,
                'n_groups': 2,
                'bcv': 0.32,
                'bcv_category': 'moderate'
            },
            warnings=[],
            all_scores={
                'DESeq2': 88.0,
                'edgeR': 82.0,
                'limma-voom': 75.0,
                'Wilcoxon': 55.0,
                'edgeR_robust': 70.0
            }
        )
        
        display_recommendation(recommendation)
        
        # Demo results
        results = {
            'timestamp': datetime.now().isoformat(),
            'raptor_version': '2.2.0',
            'module': 'M4',
            'mode': 'demo',
            'recommendation': {
                'primary_pipeline': 'DESeq2',
                'confidence': 88.0,
                'reason': 'DESeq2: appropriate for small samples, suitable for moderate dispersion',
                'alternative': 'edgeR',
                'alternative_confidence': 82.0
            },
            'decision_factors': {
                'sample_size_per_group': 6,
                'bcv': 0.32,
                'bcv_category': 'moderate'
            }
        }
        
    else:
        # Real recommendation with RAPTOR
        recommender = PipelineRecommender(profile)
        recommendation = recommender.get_recommendation()
        
        display_recommendation(recommendation)
        
        if verbose and hasattr(recommendation, 'r_code_primary'):
            print("\n" + "="*70)
            print("  📝 R CODE FOR PRIMARY RECOMMENDATION:")
            print("="*70)
            print(recommendation.r_code_primary)
            print("="*70)
        
        # Real results
        results = {
            'timestamp': datetime.now().isoformat(),
            'raptor_version': '2.2.0',
            'module': 'M4',
            'recommendation': {
                'primary_pipeline': recommendation.primary_pipeline,
                'confidence': recommendation.primary_score,
                'reason': recommendation.primary_reason,
                'alternative': recommendation.alternative_pipeline,
                'alternative_confidence': recommendation.alternative_score,
                'alternative_reason': recommendation.alternative_reason
            },
            'decision_factors': recommendation.decision_factors,
            'warnings': recommendation.warnings,
            'all_scores': recommendation.all_scores
        }
    
    # Save recommendation
    rec_file = output_path / 'recommendation.json'
    with open(rec_file, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    print(f"\n  ✓ Recommendation saved: {rec_file}")
    
    return results, output_path, recommendation


def main():
    parser = argparse.ArgumentParser(
        description='🦖 RAPTOR v2.2.0 Pipeline Recommender (Module 4)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Demo mode
  python 04_recommender.py --demo
  
  # From existing profile (recommended)
  python 04_recommender.py --profile results/profile.json
  
  # Auto-profile from counts
  python 04_recommender.py --counts results/qc/counts_clean.csv
  
  # With metadata for better profiling
  python 04_recommender.py \\
      --counts results/qc/counts_clean.csv \\
      --metadata results/quick_counts/sample_info.csv
  
  # Verbose mode (show R code)
  python 04_recommender.py --profile results/profile.json --verbose

CLI Equivalent:
  raptor recommend --profile results/profile.json

Workflow:
  Module 1: raptor quick-count → quick_gene_counts.csv
  Module 2: raptor qc → counts_clean.csv
  Module 3: raptor profile → profile.json
  Module 4: THIS SCRIPT → recommendation.json
  
  Next: Run the recommended DE pipeline in R

Supported Pipelines:
  • DESeq2      - General purpose, handles batch effects
  • edgeR       - Fast, handles low counts and high dispersion
  • limma-voom  - Very fast, best for large samples (n≥20)
  • Wilcoxon    - Non-parametric, best for very large samples (n≥8)
  • edgeR_robust - Handles outliers well

Output:
  results/recommendation.json  (Pipeline recommendation with reasoning)
        """
    )
    
    parser.add_argument('--profile', '-p',
                       help='Pre-generated profile JSON file from Module 3')
    parser.add_argument('--counts', '-c',
                       help='Count matrix CSV for auto-profiling')
    parser.add_argument('--metadata', '-m',
                       help='Sample metadata CSV (for auto-profiling)')
    parser.add_argument('--output', '-o',
                       default=DEFAULT_OUTPUT_DIR,
                       help=f'Output directory (default: {DEFAULT_OUTPUT_DIR})')
    parser.add_argument('--verbose', '-v',
                       action='store_true',
                       help='Show detailed recommendation reasoning and R code')
    parser.add_argument('--demo', action='store_true',
                       help='Run in demo mode with simulated data')
    
    args = parser.parse_args()
    
    print_banner()
    
    # Validate inputs
    if not args.demo and not args.profile and not args.counts:
        print("❌ ERROR: Either --profile, --counts, or --demo is required")
        parser.print_help()
        sys.exit(1)
    
    # Run recommendation
    results, output_dir, recommendation = run_recommendation(
        profile_file=args.profile,
        counts_file=args.counts,
        metadata_file=args.metadata,
        output_dir=args.output,
        demo=args.demo,
        verbose=args.verbose
    )
    
    # Final summary
    print("\n" + "="*70)
    print("  ✅ MODULE 4 (PIPELINE RECOMMENDATION) COMPLETE!")
    print("="*70)
    
    print(f"\n  📂 Output Directory: {output_dir}")
    print(f"\n  📊 Output Files:")
    print(f"     • recommendation.json  - Complete recommendation with reasoning")
    
    # Get recommended pipeline
    rec_pipeline = results.get('recommendation', {}).get('primary_pipeline', 'DESeq2')
    
    print(f"\n  🔜 Next Steps:")
    print(f"\n     Run Differential Expression Analysis in R:")
    print(f"     Recommended Pipeline: {rec_pipeline}")
    print(f"")
    
    # Show R code snippet for recommended pipeline
    if rec_pipeline == 'DESeq2':
        print(f"     # DESeq2 Analysis")
        print(f"     library(DESeq2)")
        print(f"     dds <- DESeqDataSetFromMatrix(counts, colData, design = ~ condition)")
        print(f"     dds <- DESeq(dds)")
        print(f"     results <- results(dds)")
    elif rec_pipeline == 'edgeR':
        print(f"     # edgeR Analysis")
        print(f"     library(edgeR)")
        print(f"     dge <- DGEList(counts)")
        print(f"     dge <- calcNormFactors(dge)")
        print(f"     dge <- estimateDisp(dge, design)")
        print(f"     fit <- glmQLFit(dge, design)")
    elif rec_pipeline == 'limma-voom':
        print(f"     # limma-voom Analysis")
        print(f"     library(limma)")
        print(f"     dge <- DGEList(counts)")
        print(f"     dge <- calcNormFactors(dge)")
        print(f"     v <- voom(dge, design)")
        print(f"     fit <- lmFit(v, design)")
        print(f"     fit <- eBayes(fit)")
    elif rec_pipeline == 'Wilcoxon':
        print(f"     # Wilcoxon Test")
        print(f"     library(edgeR)")
        print(f"     dge <- DGEList(counts)")
        print(f"     dge <- calcNormFactors(dge)")
        print(f"     # Use TMM-normalized counts with wilcox.test()")
    
    print(f"\n     For complete R code, use --verbose flag")
    
    print("\n" + "="*70)
    print("  Making free science for everybody around the world 🌍")
    print("="*70 + "\n")


if __name__ == '__main__':
    main()
