#!/usr/bin/env python3

"""
Generate QC summary report for quantification.

Creates HTML report with quantification statistics and basic QC metrics.
Supports both Salmon and Kallisto output formats.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
"""

import argparse
import os
import sys
import json
from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime
import pandas as pd
import numpy as np


# =============================================================================
# CONSTANTS
# =============================================================================

DEFAULT_OUTPUT_DIR = "results/quick_counts"
OUTPUT_COUNTS_FILE = "quick_gene_counts.csv"


def load_salmon_meta(quant_dir: str) -> Dict:
    """Load Salmon meta_info.json."""
    meta_file = Path(quant_dir) / 'aux_info' / 'meta_info.json'
    
    if not meta_file.exists():
        return {}
    
    with open(meta_file, 'r') as f:
        return json.load(f)


def load_kallisto_run_info(quant_dir: str) -> Dict:
    """Load Kallisto run_info.json."""
    run_file = Path(quant_dir) / 'run_info.json'
    
    if not run_file.exists():
        return {}
    
    with open(run_file, 'r') as f:
        return json.load(f)


def collect_quant_stats(quant_base_dir: str, tool: str = 'kallisto') -> pd.DataFrame:
    """
    Collect quantification statistics from all samples.
    
    Parameters
    ----------
    quant_base_dir : str
        Base directory containing sample folders
    tool : str
        'kallisto' or 'salmon'
    
    Returns
    -------
    pd.DataFrame
        Statistics for each sample
    """
    quant_path = Path(quant_base_dir)
    stats = []
    
    for sample_dir in sorted(quant_path.iterdir()):
        if not sample_dir.is_dir():
            continue
        
        sample_id = sample_dir.name
        
        if tool == 'kallisto':
            run_info = load_kallisto_run_info(str(sample_dir))
            
            # Load abundance to get total counts
            abundance_file = sample_dir / 'abundance.tsv'
            total_counts = 0
            n_transcripts = 0
            if abundance_file.exists():
                abundance = pd.read_csv(abundance_file, sep='\t')
                total_counts = abundance['est_counts'].sum()
                n_transcripts = len(abundance)
            
            sample_stats = {
                'sample_id': sample_id,
                'n_processed': run_info.get('n_processed', 0),
                'n_mapped': run_info.get('n_pseudoaligned', 0),
                'mapping_rate': run_info.get('p_pseudoaligned', 0),
                'n_unique': run_info.get('n_unique', 0),
                'tool_version': run_info.get('kallisto_version', 'unknown'),
                'total_counts': total_counts,
                'n_transcripts': n_transcripts,
            }
        else:  # salmon
            meta_info = load_salmon_meta(str(sample_dir))
            
            # Load quant.sf to get total counts
            quant_file = sample_dir / 'quant.sf'
            total_counts = 0
            n_transcripts = 0
            if quant_file.exists():
                quant = pd.read_csv(quant_file, sep='\t')
                total_counts = quant['NumReads'].sum()
                n_transcripts = len(quant)
            
            sample_stats = {
                'sample_id': sample_id,
                'n_processed': meta_info.get('num_processed', 0),
                'n_mapped': meta_info.get('num_mapped', 0),
                'mapping_rate': meta_info.get('percent_mapped', 0),
                'n_unique': meta_info.get('num_decoy_fragments', 0),
                'tool_version': meta_info.get('salmon_version', 'unknown'),
                'total_counts': total_counts,
                'n_transcripts': n_transcripts,
            }
        
        stats.append(sample_stats)
    
    return pd.DataFrame(stats)


def generate_html_report(
    stats_df: pd.DataFrame,
    counts_df: pd.DataFrame,
    output_path: str,
    tool: str = 'kallisto'
) -> None:
    """
    Generate HTML QC report.
    
    Parameters
    ----------
    stats_df : pd.DataFrame
        Sample statistics
    counts_df : pd.DataFrame
        Count matrix
    output_path : str
        Output HTML file path
    tool : str
        'kallisto' or 'salmon'
    """
    # Calculate metrics
    n_samples = len(stats_df)
    mean_mapping_rate = stats_df['mapping_rate'].mean() if 'mapping_rate' in stats_df.columns else 0
    total_reads = stats_df['n_processed'].sum() if 'n_processed' in stats_df.columns else 0
    
    # Count matrix stats
    n_genes = len(counts_df)
    genes_detected = (counts_df.sum(axis=1) > 0).sum()
    mean_counts = counts_df.sum().mean()
    
    # Library size stats
    library_sizes = counts_df.sum()
    lib_cv = library_sizes.std() / library_sizes.mean() * 100 if library_sizes.mean() > 0 else 0
    
    # Tool-specific naming
    tool_name = tool.capitalize()
    mapping_term = "Pseudoalignment" if tool == 'kallisto' else "Mapping"
    
    # Determine color scheme based on tool
    primary_color = "#9b59b6" if tool == 'kallisto' else "#3498db"
    gradient_start = "#667eea" if tool == 'kallisto' else "#11998e"
    gradient_end = "#764ba2" if tool == 'kallisto' else "#38ef7d"
    
    # HTML template
    html = f"""
<!DOCTYPE html>
<html>
<head>
    <title>RAPTOR Quick-Count QC Report ({tool_name})</title>
    <style>
        body {{
            font-family: 'Segoe UI', Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background: #f5f5f5;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid {primary_color};
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            margin-top: 30px;
        }}
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        .summary-card {{
            background: linear-gradient(135deg, {gradient_start} 0%, {gradient_end} 100%);
            color: white;
            padding: 20px;
            border-radius: 10px;
            text-align: center;
        }}
        .summary-card h3 {{
            margin: 0;
            font-size: 14px;
            opacity: 0.9;
        }}
        .summary-card .value {{
            font-size: 32px;
            font-weight: bold;
            margin: 10px 0;
        }}
        .summary-card.good {{
            background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%);
        }}
        .summary-card.warning {{
            background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}
        th, td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background: {primary_color};
            color: white;
        }}
        tr:hover {{
            background: #f5f5f5;
        }}
        .status-good {{
            color: #27ae60;
            font-weight: bold;
        }}
        .status-warning {{
            color: #f39c12;
            font-weight: bold;
        }}
        .status-bad {{
            color: #e74c3c;
            font-weight: bold;
        }}
        .next-steps {{
            background: #ecf0f1;
            padding: 20px;
            border-radius: 10px;
            margin-top: 20px;
        }}
        .next-steps pre {{
            background: #2c3e50;
            color: #ecf0f1;
            padding: 15px;
            border-radius: 5px;
            overflow-x: auto;
        }}
        .footer {{
            margin-top: 30px;
            padding-top: 20px;
            border-top: 1px solid #ddd;
            color: #666;
            font-size: 12px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🦖 RAPTOR Quick-Count QC Report ({tool_name})</h1>
        <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        <p>Module 1: Quantify → Output: <code>{OUTPUT_COUNTS_FILE}</code></p>
        
        <h2>📊 Summary</h2>
        <div class="summary-grid">
            <div class="summary-card">
                <h3>Samples</h3>
                <div class="value">{n_samples}</div>
            </div>
            <div class="summary-card {'good' if mean_mapping_rate > 70 else 'warning'}">
                <h3>Mean {mapping_term} Rate</h3>
                <div class="value">{mean_mapping_rate:.1f}%</div>
            </div>
            <div class="summary-card">
                <h3>Total Reads</h3>
                <div class="value">{total_reads/1e6:.1f}M</div>
            </div>
            <div class="summary-card {'good' if lib_cv < 30 else 'warning'}">
                <h3>Library Size CV</h3>
                <div class="value">{lib_cv:.1f}%</div>
            </div>
        </div>
        
        <div class="summary-grid">
            <div class="summary-card">
                <h3>Genes in Matrix</h3>
                <div class="value">{n_genes:,}</div>
            </div>
            <div class="summary-card good">
                <h3>Genes Detected</h3>
                <div class="value">{genes_detected:,}</div>
            </div>
            <div class="summary-card">
                <h3>Mean Library Size</h3>
                <div class="value">{mean_counts/1e6:.1f}M</div>
            </div>
        </div>
        
        <h2>📋 Per-Sample Statistics</h2>
        <table>
            <thead>
                <tr>
                    <th>Sample</th>
                    <th>Processed Reads</th>
                    <th>Mapped</th>
                    <th>{mapping_term} Rate</th>
                    <th>Library Size</th>
                    <th>Status</th>
                </tr>
            </thead>
            <tbody>
"""
    
    # Add sample rows
    for _, row in stats_df.iterrows():
        sample_id = row['sample_id']
        lib_size = library_sizes.get(sample_id, 0) if sample_id in library_sizes.index else 0
        mapping_rate = row.get('mapping_rate', 0)
        
        # Determine status
        if mapping_rate >= 70 and lib_size >= 1e6:
            status_class = 'status-good'
            status = '✓ Good'
        elif mapping_rate >= 50 and lib_size >= 0.5e6:
            status_class = 'status-warning'
            status = '⚠ Check'
        else:
            status_class = 'status-bad'
            status = '✗ Poor'
        
        html += f"""
                <tr>
                    <td>{sample_id}</td>
                    <td>{row.get('n_processed', 0):,.0f}</td>
                    <td>{row.get('n_mapped', 0):,.0f}</td>
                    <td>{mapping_rate:.1f}%</td>
                    <td>{lib_size:,.0f}</td>
                    <td class="{status_class}">{status}</td>
                </tr>
"""
    
    html += f"""
            </tbody>
        </table>
        
        <div class="next-steps">
            <h2>📋 Next Steps (RAPTOR Workflow)</h2>
            <p>This report provides basic QC metrics from {tool_name} quantification.</p>
            <p>Continue with the RAPTOR workflow:</p>
            <pre>
# Module 2: Quality Assessment & Outlier Detection
raptor qc --counts {DEFAULT_OUTPUT_DIR}/{OUTPUT_COUNTS_FILE}

# Module 3: Data Profiling (32 features)
raptor profile --counts {DEFAULT_OUTPUT_DIR}/{OUTPUT_COUNTS_FILE}

# Module 4: Get Pipeline Recommendation
raptor recommend
            </pre>
        </div>
        
        <div class="footer">
            <p>Generated by RAPTOR v2.2.0 - RNA-seq Analysis Pipeline Testing & Optimization Resource</p>
            <p>Author: Ayeh Bolouki | Email: ayehbolouki1988@gmail.com</p>
            <p>GitHub: <a href="https://github.com/AyehBlk/RAPTOR">https://github.com/AyehBlk/RAPTOR</a></p>
        </div>
    </div>
</body>
</html>
"""
    
    # Write report
    with open(output_path, 'w') as f:
        f.write(html)
    
    print(f"✓ Report saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate QC report for quantification (RAPTOR v2.2.0)'
    )
    
    parser.add_argument(
        '--quant-dir', '-q',
        required=True,
        help='Directory containing sample quant folders'
    )
    
    parser.add_argument(
        '--counts', '-c',
        required=True,
        help='Count matrix CSV file'
    )
    
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output HTML report path'
    )
    
    parser.add_argument(
        '--tool',
        choices=['kallisto', 'salmon'],
        default='kallisto',
        help='Quantification tool (default: kallisto)'
    )
    
    args = parser.parse_args()
    
    print(f"📊 Collecting {args.tool} quantification statistics...")
    stats_df = collect_quant_stats(args.quant_dir, args.tool)
    
    print(f"📂 Loading count matrix...")
    counts_df = pd.read_csv(args.counts, index_col=0)
    
    print(f"📝 Generating HTML report...")
    generate_html_report(stats_df, counts_df, args.output, args.tool)
    
    print(f"\n✅ Report generated successfully!")


if __name__ == '__main__':
    main()
