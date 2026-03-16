#!/usr/bin/env python3

"""
Combine quantification outputs into gene-level count matrix.

Reads individual abundance files (Kallisto abundance.tsv or Salmon quant.sf)
and aggregates to gene level using tx2gene mapping.

Output Files (Architecture Compliant):
- results/quick_counts/quick_gene_counts.csv
- results/quick_counts/quick_tpm.csv
- results/quick_counts/effective_lengths.csv

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
"""

import argparse
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pandas as pd
import numpy as np


# =============================================================================
# CONSTANTS - Architecture Compliant
# =============================================================================

DEFAULT_OUTPUT_DIR = "results/quick_counts"
OUTPUT_COUNTS_FILE = "quick_gene_counts.csv"
OUTPUT_TPM_FILE = "quick_tpm.csv"
OUTPUT_LENGTHS_FILE = "effective_lengths.csv"


def load_tx2gene(tx2gene_path: str) -> Dict[str, str]:
    """
    Load transcript-to-gene mapping.
    
    Supports:
    - TSV/CSV with columns: transcript_id, gene_id
    - GTF file (extracts transcript_id and gene_id attributes)
    
    Parameters
    ----------
    tx2gene_path : str
        Path to tx2gene mapping file
    
    Returns
    -------
    Dict[str, str]
        Mapping from transcript_id to gene_id
    """
    path = Path(tx2gene_path)
    
    if path.suffix.lower() in ['.gtf', '.gff', '.gff3']:
        return _load_tx2gene_from_gtf(tx2gene_path)
    else:
        return _load_tx2gene_from_table(tx2gene_path)


def _load_tx2gene_from_table(path: str) -> Dict[str, str]:
    """Load tx2gene from TSV/CSV file."""
    # Try different separators
    for sep in ['\t', ',']:
        try:
            df = pd.read_csv(path, sep=sep)
            break
        except:
            continue
    
    # Normalize column names
    df.columns = [c.lower().strip() for c in df.columns]
    
    # Find transcript and gene columns
    tx_col = None
    gene_col = None
    
    tx_aliases = ['transcript_id', 'txname', 'transcript', 'tx_id', 'target_id', 'ensembl_transcript_id']
    gene_aliases = ['gene_id', 'geneid', 'gene', 'gene_name', 'ens_gene', 'ensembl_gene_id']
    
    for alias in tx_aliases:
        if alias in df.columns:
            tx_col = alias
            break
    
    for alias in gene_aliases:
        if alias in df.columns:
            gene_col = alias
            break
    
    # Fallback: assume first two columns
    if tx_col is None or gene_col is None:
        df.columns = ['transcript_id', 'gene_id'] + list(df.columns[2:])
        tx_col = 'transcript_id'
        gene_col = 'gene_id'
    
    return dict(zip(df[tx_col], df[gene_col]))


def _load_tx2gene_from_gtf(path: str) -> Dict[str, str]:
    """Extract tx2gene mapping from GTF file."""
    import re
    tx2gene = {}
    
    with open(path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            feature_type = fields[2]
            if feature_type not in ['transcript', 'exon']:
                continue
            
            attributes = fields[8]
            
            # Extract transcript_id
            tx_match = None
            for pattern in [r'transcript_id "([^"]+)"', r'transcript_id=([^;]+)']:
                match = re.search(pattern, attributes)
                if match:
                    tx_match = match.group(1)
                    break
            
            # Extract gene_id
            gene_match = None
            for pattern in [r'gene_id "([^"]+)"', r'gene_id=([^;]+)']:
                match = re.search(pattern, attributes)
                if match:
                    gene_match = match.group(1)
                    break
            
            if tx_match and gene_match:
                tx2gene[tx_match] = gene_match
    
    return tx2gene


def load_abundance(quant_dir: str, tool: str = 'kallisto') -> pd.DataFrame:
    """
    Load abundance file from Kallisto or Salmon.
    
    Parameters
    ----------
    quant_dir : str
        Directory containing abundance file
    tool : str
        'kallisto' or 'salmon'
    
    Returns
    -------
    pd.DataFrame
        Abundance data
    """
    quant_path = Path(quant_dir)
    
    if tool == 'kallisto':
        abundance_file = quant_path / 'abundance.tsv'
        if not abundance_file.exists():
            raise FileNotFoundError(f"abundance.tsv not found in: {quant_dir}")
        df = pd.read_csv(abundance_file, sep='\t')
        # Kallisto columns: target_id, length, eff_length, est_counts, tpm
        df = df.rename(columns={'target_id': 'Name', 'est_counts': 'NumReads', 'tpm': 'TPM'})
    else:  # salmon
        abundance_file = quant_path / 'quant.sf'
        if not abundance_file.exists():
            raise FileNotFoundError(f"quant.sf not found in: {quant_dir}")
        df = pd.read_csv(abundance_file, sep='\t')
        # Salmon columns: Name, Length, EffectiveLength, TPM, NumReads
    
    return df


def aggregate_to_gene(
    abundance_df: pd.DataFrame,
    tx2gene: Dict[str, str],
    count_col: str = 'NumReads'
) -> pd.Series:
    """
    Aggregate transcript counts to gene level.
    
    Parameters
    ----------
    abundance_df : pd.DataFrame
        Abundance data with 'Name' column
    tx2gene : Dict[str, str]
        Transcript to gene mapping
    count_col : str
        Column to aggregate
    
    Returns
    -------
    pd.Series
        Gene-level values
    """
    abundance_df = abundance_df.copy()
    abundance_df['gene_id'] = abundance_df['Name'].map(tx2gene)
    
    # Report unmapped
    unmapped = abundance_df['gene_id'].isna().sum()
    if unmapped > 0:
        print(f"  Warning: {unmapped} transcripts not found in tx2gene mapping")
        abundance_df = abundance_df.dropna(subset=['gene_id'])
    
    # Aggregate
    gene_values = abundance_df.groupby('gene_id')[count_col].sum()
    
    return gene_values


def combine_samples(
    quant_dir: str,
    tx2gene: Dict[str, str],
    tool: str = 'kallisto',
    sample_info_path: Optional[str] = None
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Combine all quantification outputs into count matrices.
    
    Parameters
    ----------
    quant_dir : str
        Directory containing sample subdirectories
    tx2gene : Dict[str, str]
        Transcript to gene mapping
    tool : str
        'kallisto' or 'salmon'
    sample_info_path : str, optional
        Path to sample info CSV for ordering
    
    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
        (counts_df, tpm_df, length_df)
    """
    quant_path = Path(quant_dir)
    
    # Get sample directories
    if sample_info_path:
        sample_info = pd.read_csv(sample_info_path)
        sample_ids = sample_info['sample_id'].tolist()
        sample_dirs = [quant_path / sid for sid in sample_ids]
    else:
        sample_dirs = sorted([d for d in quant_path.iterdir() if d.is_dir()])
        sample_ids = [d.name for d in sample_dirs]
    
    print(f"Combining {len(sample_ids)} samples...")
    
    counts_dict = {}
    tpm_dict = {}
    length_dict = {}
    
    for sample_id, sample_dir in zip(sample_ids, sample_dirs):
        if not sample_dir.exists():
            print(f"  Warning: Directory not found for {sample_id}")
            continue
        
        print(f"  Processing: {sample_id}")
        
        try:
            abundance_df = load_abundance(str(sample_dir), tool)
            
            counts_dict[sample_id] = aggregate_to_gene(abundance_df, tx2gene, 'NumReads')
            tpm_dict[sample_id] = aggregate_to_gene(abundance_df, tx2gene, 'TPM')
            
            # Effective length (use mean for gene level)
            if 'EffectiveLength' in abundance_df.columns:
                length_col = 'EffectiveLength'
            elif 'eff_length' in abundance_df.columns:
                length_col = 'eff_length'
            else:
                length_col = None
            
            if length_col:
                abundance_df['gene_id'] = abundance_df['Name'].map(tx2gene)
                length_dict[sample_id] = abundance_df.groupby('gene_id')[length_col].mean()
                
        except Exception as e:
            print(f"  Error processing {sample_id}: {e}")
            continue
    
    if not counts_dict:
        raise ValueError("No samples could be processed!")
    
    # Create DataFrames
    counts_df = pd.DataFrame(counts_dict).fillna(0)
    tpm_df = pd.DataFrame(tpm_dict).fillna(0)
    length_df = pd.DataFrame(length_dict).fillna(0) if length_dict else pd.DataFrame()
    
    # Round counts to integers
    counts_df = counts_df.round().astype(int)
    
    return counts_df, tpm_df, length_df


def save_matrices(
    counts_df: pd.DataFrame,
    tpm_df: pd.DataFrame,
    length_df: pd.DataFrame,
    output_dir: str
) -> None:
    """
    Save count matrices to files with architecture-compliant names.
    
    Parameters
    ----------
    counts_df : pd.DataFrame
        Gene count matrix
    tpm_df : pd.DataFrame
        TPM matrix
    length_df : pd.DataFrame
        Effective length matrix
    output_dir : str
        Output directory
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Save counts with architecture-compliant name
    counts_file = output_path / OUTPUT_COUNTS_FILE
    counts_df.to_csv(counts_file)
    print(f"✓ Saved: {counts_file}")
    
    # Also save as TSV for compatibility
    counts_df.to_csv(output_path / 'quick_gene_counts.tsv', sep='\t')
    
    # Save TPM
    tpm_file = output_path / OUTPUT_TPM_FILE
    tpm_df.to_csv(tpm_file)
    print(f"✓ Saved: {tpm_file}")
    
    # Save effective lengths
    if not length_df.empty:
        length_file = output_path / OUTPUT_LENGTHS_FILE
        length_df.to_csv(length_file)
        print(f"✓ Saved: {length_file}")
    
    # Summary statistics
    summary = {
        'n_genes': len(counts_df),
        'n_samples': len(counts_df.columns),
        'total_counts': int(counts_df.sum().sum()),
        'mean_counts_per_sample': int(counts_df.sum().mean()),
        'genes_detected': int((counts_df.sum(axis=1) > 0).sum()),
    }
    
    summary_file = output_path / 'count_summary.txt'
    with open(summary_file, 'w') as f:
        f.write("RAPTOR Quick-Count Summary\n")
        f.write("=" * 40 + "\n")
        f.write(f"Output: {OUTPUT_COUNTS_FILE}\n")
        f.write("=" * 40 + "\n")
        for key, value in summary.items():
            f.write(f"{key}: {value:,}\n")
    
    print(f"✓ Saved: {summary_file}")
    
    # Print summary
    print(f"\n📊 Count Matrix Summary:")
    print(f"   Genes: {summary['n_genes']:,}")
    print(f"   Samples: {summary['n_samples']}")
    print(f"   Total counts: {summary['total_counts']:,}")
    print(f"   Genes detected: {summary['genes_detected']:,}")


def main():
    parser = argparse.ArgumentParser(
        description='Combine quantification outputs into gene-level count matrix (RAPTOR v2.2.0)'
    )
    
    parser.add_argument(
        '--quant-dir', '-q',
        required=True,
        help='Directory containing sample quant folders'
    )
    
    parser.add_argument(
        '--tx2gene', '-t',
        required=True,
        help='Transcript-to-gene mapping file (TSV/CSV/GTF)'
    )
    
    parser.add_argument(
        '--output', '-o',
        default=DEFAULT_OUTPUT_DIR,
        help=f'Output directory (default: {DEFAULT_OUTPUT_DIR})'
    )
    
    parser.add_argument(
        '--sample-info', '-s',
        help='Sample info CSV (optional, for consistent ordering)'
    )
    
    parser.add_argument(
        '--tool',
        choices=['kallisto', 'salmon'],
        default='kallisto',
        help='Quantification tool used (default: kallisto)'
    )
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.isdir(args.quant_dir):
        print(f"❌ Error: Quant directory not found: {args.quant_dir}")
        sys.exit(1)
    
    if not os.path.isfile(args.tx2gene):
        print(f"❌ Error: tx2gene file not found: {args.tx2gene}")
        sys.exit(1)
    
    # Load tx2gene
    print(f"📂 Loading tx2gene mapping: {args.tx2gene}")
    tx2gene = load_tx2gene(args.tx2gene)
    print(f"   Loaded {len(tx2gene):,} transcript-gene mappings")
    
    # Combine samples
    counts_df, tpm_df, length_df = combine_samples(
        args.quant_dir,
        tx2gene,
        args.tool,
        args.sample_info
    )
    
    # Save matrices
    save_matrices(counts_df, tpm_df, length_df, args.output)
    
    print(f"\n✅ Count matrices generated successfully!")
    print(f"\n📋 Next step (Module 2):")
    print(f"   raptor qc --counts {args.output}/{OUTPUT_COUNTS_FILE}")


if __name__ == '__main__':
    main()
