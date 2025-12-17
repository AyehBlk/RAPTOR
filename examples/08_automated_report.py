#!/usr/bin/env python3
"""
RAPTOR v2.1.1 Example Script: Automated Reporting

Demonstrates publication-ready report generation with:
- Executive summary
- Data quality assessment
- Statistical results interpretation
- Biological context
- Comprehensive visualizations
- Actionable recommendations
- Multiple output formats (HTML, Markdown, PDF)

Author: Ayeh Bolouki
License: MIT
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

# Check for dependencies
try:
    import numpy as np
    import pandas as pd
except ImportError:
    print("ERROR: numpy and pandas are required")
    print("Install with: pip install numpy pandas")
    sys.exit(1)

# RAPTOR imports with fallback
RAPTOR_AVAILABLE = True
try:
    from raptor.automated_reporting import AutomatedReporter, generate_analysis_report
except ImportError:
    RAPTOR_AVAILABLE = False
    print("NOTE: RAPTOR modules not available. Running in demo mode.")


def print_banner():
    """Print RAPTOR banner."""
    print("""
    ╔══════════════════════════════════════════════════════════════╗
    ║       🦖 RAPTOR v2.1.1 - Automated Reporting                ║
    ║                                                              ║
    ║   Publication-Ready Reports with Biological Interpretation  ║
    ╚══════════════════════════════════════════════════════════════╝
    """)


def generate_demo_results():
    """Generate demonstration analysis results."""
    np.random.seed(42)
    
    # Generate DE genes
    n_de = 150
    de_genes = []
    for i in range(n_de):
        lfc = np.random.normal(0, 2)
        pval = np.random.uniform(0.00001, 0.05)
        de_genes.append({
            'gene_id': f'Gene_{i+1:05d}',
            'gene_symbol': f'GENE{i+1}',
            'log2fc': lfc,
            'pvalue': pval,
            'padj': pval * 1.5,
            'direction': 'up' if lfc > 0 else 'down'
        })
    
    return {
        'de_genes': de_genes,
        'n_total_genes': 20000,
        'n_de_genes': n_de,
        'n_upregulated': sum(1 for g in de_genes if g['direction'] == 'up'),
        'n_downregulated': sum(1 for g in de_genes if g['direction'] == 'down'),
        'top_pathways': [
            {'name': 'Cell cycle regulation', 'pvalue': 0.0001, 'genes': 15},
            {'name': 'Apoptosis', 'pvalue': 0.0005, 'genes': 12},
            {'name': 'DNA repair', 'pvalue': 0.001, 'genes': 10},
            {'name': 'Metabolism', 'pvalue': 0.005, 'genes': 8}
        ],
        'pipeline': 'Salmon-edgeR',
        'pipeline_id': 3,
        'analysis_date': datetime.now().isoformat()
    }


def generate_demo_profile():
    """Generate demonstration data profile."""
    return {
        'design': {
            'n_samples': 12,
            'n_genes': 20000,
            'n_conditions': 2,
            'samples_per_condition': 6,
            'is_paired': False
        },
        'library_stats': {
            'mean': 25000000.0,
            'cv': 0.18
        },
        'count_distribution': {
            'zero_pct': 42.5
        },
        'biological_variation': {
            'bcv': 0.35
        },
        'sequencing': {
            'depth_category': 'high',
            'mean_depth': 25000000.0
        },
        'complexity': {
            'score': 65.0
        },
        'quality_metrics': {
            'overall_score': 78.5,
            'status': 'good',
            'components': {
                'library_quality': {'score': 85, 'status': 'good'},
                'gene_detection': {'score': 75, 'status': 'good'},
                'batch_effects': {'score': 70, 'status': 'warning', 'flags': ['Minor batch effect detected']}
            }
        },
        'batch_effects': {
            'batch_detected': True,
            'batch_variable': 'batch'
        }
    }


def generate_html_report(results, profile, output_dir):
    """Generate HTML report in demo mode."""
    
    n_de = results['n_de_genes']
    n_up = results['n_upregulated']
    n_down = results['n_downregulated']
    
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>🦖 RAPTOR Analysis Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; background: #f5f5f5; }}
        .container {{ max-width: 1200px; margin: 0 auto; background: white; padding: 40px; border-radius: 10px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }}
        h1 {{ color: #2E7D32; border-bottom: 3px solid #2E7D32; padding-bottom: 10px; }}
        h2 {{ color: #1976D2; margin-top: 30px; border-bottom: 1px solid #ddd; padding-bottom: 5px; }}
        .summary-box {{ background: linear-gradient(135deg, #E8F5E9 0%, #C8E6C9 100%); padding: 25px; border-radius: 10px; margin: 20px 0; border-left: 5px solid #2E7D32; }}
        .warning-box {{ background: linear-gradient(135deg, #FFF3E0 0%, #FFE0B2 100%); padding: 20px; border-radius: 10px; margin: 20px 0; border-left: 5px solid #FF9800; }}
        .stat-grid {{ display: grid; grid-template-columns: repeat(4, 1fr); gap: 20px; margin: 20px 0; }}
        .stat-box {{ background: #f8f9fa; padding: 20px; border-radius: 8px; text-align: center; }}
        .stat-box h3 {{ color: #1976D2; margin: 0; font-size: 2em; }}
        .stat-box p {{ color: #666; margin: 5px 0 0 0; }}
        table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
        th {{ background-color: #2E7D32; color: white; }}
        tr:nth-child(even) {{ background-color: #f8f9fa; }}
        .badge {{ display: inline-block; padding: 4px 12px; border-radius: 15px; font-size: 0.85em; }}
        .badge-up {{ background: #C8E6C9; color: #2E7D32; }}
        .badge-down {{ background: #FFCDD2; color: #C62828; }}
        .footer {{ margin-top: 40px; padding-top: 20px; border-top: 1px solid #ddd; text-align: center; color: #666; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🦖 RAPTOR RNA-seq Analysis Report</h1>
        <p><strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
        <p><strong>Pipeline:</strong> {results['pipeline']} (Pipeline {results['pipeline_id']})</p>
        
        <div class="summary-box">
            <h2 style="margin-top: 0; border: none;">Executive Summary</h2>
            <ul>
                <li>Identified <strong>{n_de} differentially expressed genes</strong></li>
                <li><strong>{n_up} upregulated</strong>, <strong>{n_down} downregulated</strong></li>
                <li>Data quality score: <strong>{profile['quality_metrics']['overall_score']:.1f}/100</strong></li>
                <li>Analyzed {profile['design']['n_samples']} samples across {profile['design']['n_conditions']} conditions</li>
            </ul>
        </div>
        
        <div class="stat-grid">
            <div class="stat-box">
                <h3>{n_de}</h3>
                <p>DE Genes</p>
            </div>
            <div class="stat-box">
                <h3>{n_up}</h3>
                <p>Upregulated</p>
            </div>
            <div class="stat-box">
                <h3>{n_down}</h3>
                <p>Downregulated</p>
            </div>
            <div class="stat-box">
                <h3>{profile['quality_metrics']['overall_score']:.0f}</h3>
                <p>Quality Score</p>
            </div>
        </div>
        
        <h2>📊 Data Overview</h2>
        <table>
            <tr><th>Metric</th><th>Value</th></tr>
            <tr><td>Total Samples</td><td>{profile['design']['n_samples']}</td></tr>
            <tr><td>Total Genes</td><td>{profile['design']['n_genes']:,}</td></tr>
            <tr><td>Mean Library Size</td><td>{profile['library_stats']['mean']/1e6:.1f}M reads</td></tr>
            <tr><td>Biological Variation (BCV)</td><td>{profile['biological_variation']['bcv']:.3f}</td></tr>
            <tr><td>Sequencing Depth</td><td>{profile['sequencing']['depth_category'].title()}</td></tr>
        </table>
        
        <h2>🔬 Top Differentially Expressed Genes</h2>
        <table>
            <tr><th>Gene ID</th><th>Symbol</th><th>Log2 FC</th><th>Adjusted P-value</th><th>Direction</th></tr>
"""
    
    # Add top 10 genes
    for gene in sorted(results['de_genes'], key=lambda x: abs(x['log2fc']), reverse=True)[:10]:
        direction = 'up' if gene['log2fc'] > 0 else 'down'
        badge_class = 'badge-up' if direction == 'up' else 'badge-down'
        html_content += f"""            <tr>
                <td>{gene['gene_id']}</td>
                <td>{gene['gene_symbol']}</td>
                <td>{gene['log2fc']:.2f}</td>
                <td>{gene['padj']:.2e}</td>
                <td><span class="badge {badge_class}">{'↑ Up' if direction == 'up' else '↓ Down'}</span></td>
            </tr>
"""
    
    html_content += """        </table>
        
        <h2>🧬 Pathway Enrichment</h2>
        <table>
            <tr><th>Pathway</th><th>P-value</th><th>Genes</th></tr>
"""
    
    for pathway in results['top_pathways']:
        html_content += f"""            <tr>
                <td>{pathway['name']}</td>
                <td>{pathway['pvalue']:.4f}</td>
                <td>{pathway['genes']}</td>
            </tr>
"""
    
    html_content += f"""        </table>
        
        <h2>📋 Methods</h2>
        <p>
            RNA-seq differential expression analysis was performed using the {results['pipeline']} pipeline.
            Raw counts were filtered to remove lowly expressed genes (counts < 10 in fewer than 
            {profile['design']['samples_per_condition']} samples). Library size normalization was applied using 
            TMM (trimmed mean of M-values). Differential expression testing was performed with 
            an FDR threshold of 0.05 and no log2 fold-change threshold.
        </p>
        <p>
            Analysis was conducted using RAPTOR v2.1.1 (RNA-seq Analysis Pipeline Testing and 
            Optimization Resource). Quality assessment indicated {profile['quality_metrics']['status']} 
            data quality with an overall score of {profile['quality_metrics']['overall_score']:.1f}/100.
        </p>
        
        <div class="footer">
            <p>🦖 Generated by RAPTOR v2.1.1</p>
            <p>Making free science for everybody around the world 🌍</p>
            <p><a href="https://github.com/AyehBlk/RAPTOR">https://github.com/AyehBlk/RAPTOR</a></p>
        </div>
    </div>
</body>
</html>
"""
    
    output_file = Path(output_dir) / 'analysis_report.html'
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    return output_file


def display_report_summary(results, profile, format_type):
    """Display report generation summary."""
    
    print("\n" + "="*70)
    print("  🦖 AUTOMATED REPORT GENERATION")
    print("="*70)
    
    # Key findings
    print("\n  📊 Key Findings:")
    print(f"     DE Genes:       {results['n_de_genes']}")
    print(f"     Upregulated:    {results['n_upregulated']}")
    print(f"     Downregulated:  {results['n_downregulated']}")
    print(f"     Quality Score:  {profile['quality_metrics']['overall_score']:.1f}/100")
    
    # Top pathways
    if 'top_pathways' in results:
        print("\n  🧬 Top Enriched Pathways:")
        for pathway in results['top_pathways'][:3]:
            print(f"     • {pathway['name']} (p={pathway['pvalue']:.4f}, {pathway['genes']} genes)")
    
    # Top DE genes
    print("\n  🔝 Top DE Genes (by fold change):")
    top_genes = sorted(results['de_genes'], key=lambda x: abs(x['log2fc']), reverse=True)[:5]
    for gene in top_genes:
        direction = '↑' if gene['log2fc'] > 0 else '↓'
        print(f"     {direction} {gene['gene_symbol']:12s} LFC: {gene['log2fc']:>6.2f}")


def run_automated_report(results_file=None, profile_file=None, format_type='html',
                         title=None, pathway=False, output_dir='report', demo=False):
    """Run automated report generation."""
    
    if demo or not RAPTOR_AVAILABLE:
        print("\n🎮 Running in DEMO mode with simulated results...")
        results = generate_demo_results()
        profile = generate_demo_profile()
        
        print(f"\n📊 Demo results generated:")
        print(f"   DE Genes: {results['n_de_genes']}")
        print(f"   Upregulated: {results['n_upregulated']}")
        print(f"   Downregulated: {results['n_downregulated']}")
        
        # Generate report
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        if format_type in ['html', 'all']:
            html_file = generate_html_report(results, profile, output_dir)
            print(f"\n📄 HTML report generated: {html_file}")
        
    else:
        # Load results
        print(f"\n📂 Loading results from: {results_file}")
        with open(results_file, 'r') as f:
            results = json.load(f)
        
        # Load profile
        print(f"📂 Loading profile from: {profile_file}")
        with open(profile_file, 'r') as f:
            profile = json.load(f)
        
        # Generate report using RAPTOR
        print(f"\n📝 Generating {format_type} report...")
        
        reporter = AutomatedReporter(results, profile, output_dir)
        report_files = reporter.generate_complete_report(format=format_type)
        
        print(f"\n📄 Reports generated:")
        for file_type, file_path in report_files.items():
            print(f"   {file_type}: {file_path}")
    
    # Display summary
    display_report_summary(results, profile, format_type)
    
    # Prepare output summary
    output = {
        'timestamp': datetime.now().isoformat(),
        'raptor_version': '2.1.1',
        'format': format_type,
        'output_dir': output_dir,
        'summary': {
            'n_de_genes': results['n_de_genes'],
            'n_upregulated': results['n_upregulated'],
            'n_downregulated': results['n_downregulated'],
            'quality_score': profile['quality_metrics']['overall_score']
        }
    }
    
    return output


def main():
    parser = argparse.ArgumentParser(
        description='🦖 RAPTOR Automated Report Generation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate HTML report from results
  python 08_automated_report.py --results de_results.json --profile data_profile.json
  
  # Generate multiple formats
  python 08_automated_report.py --results results.json --profile profile.json --format all
  
  # Custom title and output directory
  python 08_automated_report.py --results results.json --profile profile.json \\
      --title "My RNA-seq Analysis" --output my_report/
  
  # Demo mode (no data required)
  python 08_automated_report.py --demo

Output Formats:
  html     - Interactive HTML report (default)
  markdown - Markdown for documentation
  pdf      - PDF for publication
  all      - Generate all formats
        """
    )
    
    parser.add_argument('--results', help='Analysis results JSON file')
    parser.add_argument('--profile', help='Data profile JSON file')
    parser.add_argument('--format', choices=['html', 'markdown', 'pdf', 'all'],
                        default='html', help='Output format')
    parser.add_argument('--title', help='Report title')
    parser.add_argument('--pathway', action='store_true', help='Include pathway enrichment')
    parser.add_argument('--output', '-o', default='report', help='Output directory')
    parser.add_argument('--demo', action='store_true', help='Run in demo mode')
    
    args = parser.parse_args()
    
    print_banner()
    
    # Validate inputs
    if not args.demo and (not args.results or not args.profile):
        print("ERROR: --results and --profile are required, or use --demo")
        parser.print_help()
        sys.exit(1)
    
    if args.results and not Path(args.results).exists():
        print(f"ERROR: Results file not found: {args.results}")
        sys.exit(1)
    
    # Run report generation
    output = run_automated_report(
        results_file=args.results,
        profile_file=args.profile,
        format_type=args.format,
        title=args.title,
        pathway=args.pathway,
        output_dir=args.output,
        demo=args.demo
    )
    
    # Save summary
    summary_file = Path(args.output) / 'report_summary.json'
    with open(summary_file, 'w') as f:
        json.dump(output, f, indent=2, default=str)
    
    print(f"\n💾 Summary saved to: {summary_file}")
    print("\n" + "="*70)
    print("  Making free science for everybody around the world 🌍")
    print("="*70 + "\n")


if __name__ == '__main__':
    main()
