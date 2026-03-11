#!/usr/bin/env python3
"""
RAPTOR v2.2.0 - Pipeline Registry & Metadata
==============================================

Central registry for all RNA-seq quantification pipelines with comprehensive
metadata, lazy loading, and helper functions.

✅ PRODUCTION PIPELINES (Module 5) - 6/6 Complete:
===================================================
1. hisat2_featurecounts - HISAT2 + featureCounts (16GB, BAM files)
2. kallisto             - Kallisto pseudo-alignment (4GB, fastest)
3. salmon               - Salmon pseudo-alignment (8GB, best balance) ⭐
4. star_featurecounts   - STAR + featureCounts (32GB, gold standard genes)
5. star_rsem            - STAR + RSEM (32GB, gold standard isoforms)
6. star_salmon          - STAR + Salmon (32GB, BAM + bootstraps) 🌟

⏳ QUICK PIPELINES (Module 1) - Optional:
=========================================
- quick_kallisto: Fast Kallisto for QC profiling
- quick_salmon: Fast Salmon for QC profiling

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
Date: January 2026
"""

from typing import Dict, List, Optional, Type, Union
from pathlib import Path

__version__ = "2.2.0"

# =============================================================================
# BASE CLASS IMPORTS
# =============================================================================

from .base import (
    BasePipeline,
    PipelineConfig,
    PipelineResult,
    SampleInfo,
    ToolDependency,
    DependencyHandler,
    auto_detect_samples,
    create_sample_sheet
)

# =============================================================================
# PIPELINE REGISTRY (Lazy Loading)
# =============================================================================

_PIPELINE_REGISTRY: Dict[str, Type[BasePipeline]] = {}
_PIPELINE_PARAMS: Dict[str, Dict] = {}


def _load_pipelines():
    """
    Lazy load pipeline classes.
    
    This function loads all production pipelines with proper error handling.
    Pipelines are loaded on first access to minimize import time.
    """
    global _PIPELINE_REGISTRY, _PIPELINE_PARAMS
    
    if _PIPELINE_REGISTRY:
        return
    
    # =========================================================================
    # PRODUCTION PIPELINES - ALL IMPLEMENTED v2.2.0
    # =========================================================================
    
    # 1. HISAT2 + featureCounts - Low memory alternative
    try:
        from .hisat2_featurecounts.pipeline import Hisat2FeatureCountsPipeline, HISAT2_FC_CLI_PARAMS
        _PIPELINE_REGISTRY['hisat2_featurecounts'] = Hisat2FeatureCountsPipeline
        _PIPELINE_PARAMS['hisat2_featurecounts'] = HISAT2_FC_CLI_PARAMS
    except ImportError as e:
        import warnings
        warnings.warn(f"Could not load HISAT2+featureCounts pipeline: {e}")
    
    # 2. Kallisto - Fastest pseudo-alignment
    try:
        from .kallisto.pipeline import KallistoPipeline, KALLISTO_CLI_PARAMS
        _PIPELINE_REGISTRY['kallisto'] = KallistoPipeline
        _PIPELINE_PARAMS['kallisto'] = KALLISTO_CLI_PARAMS
    except ImportError as e:
        import warnings
        warnings.warn(f"Could not load Kallisto pipeline: {e}")
    
    # 3. Salmon - Best balance speed/accuracy
    try:
        from .salmon.pipeline import SalmonPipeline, SALMON_CLI_PARAMS
        _PIPELINE_REGISTRY['salmon'] = SalmonPipeline
        _PIPELINE_PARAMS['salmon'] = SALMON_CLI_PARAMS
    except ImportError as e:
        import warnings
        warnings.warn(f"Could not load Salmon pipeline: {e}")
    
    # 4. STAR + featureCounts - Gold standard for genes
    try:
        from .star_featurecounts.pipeline import StarFeatureCountsPipeline, STAR_FC_CLI_PARAMS
        _PIPELINE_REGISTRY['star_featurecounts'] = StarFeatureCountsPipeline
        _PIPELINE_PARAMS['star_featurecounts'] = STAR_FC_CLI_PARAMS
    except ImportError as e:
        import warnings
        warnings.warn(f"Could not load STAR+featureCounts pipeline: {e}")
    
    # 5. STAR + RSEM - Gold standard for isoforms
    try:
        from .star_rsem.pipeline import StarRsemPipeline, STAR_RSEM_CLI_PARAMS
        _PIPELINE_REGISTRY['star_rsem'] = StarRsemPipeline
        _PIPELINE_PARAMS['star_rsem'] = STAR_RSEM_CLI_PARAMS
    except ImportError as e:
        import warnings
        warnings.warn(f"Could not load STAR+RSEM pipeline: {e}")
    
    # 6. STAR + Salmon - Hybrid: BAM + bootstraps
    try:
        from .star_salmon.pipeline import StarSalmonPipeline, STAR_SALMON_CLI_PARAMS
        _PIPELINE_REGISTRY['star_salmon'] = StarSalmonPipeline
        _PIPELINE_PARAMS['star_salmon'] = STAR_SALMON_CLI_PARAMS
    except ImportError as e:
        import warnings
        warnings.warn(f"Could not load STAR+Salmon pipeline: {e}")
    
    # =========================================================================
    # QUICK PIPELINES (Module 1) - Optional QC/Profiling Tools
    # =========================================================================
    
    try:
        from .quick_kallisto import QuickKallistoPipeline, QUICK_KALLISTO_CLI_PARAMS
        _PIPELINE_REGISTRY['quick_kallisto'] = QuickKallistoPipeline
        _PIPELINE_PARAMS['quick_kallisto'] = QUICK_KALLISTO_CLI_PARAMS
    except ImportError:
        pass  # Optional - Quick pipelines
    
    try:
        from .quick_salmon import QuickSalmonPipeline, QUICK_SALMON_CLI_PARAMS
        _PIPELINE_REGISTRY['quick_salmon'] = QuickSalmonPipeline
        _PIPELINE_PARAMS['quick_salmon'] = QUICK_SALMON_CLI_PARAMS
    except ImportError:
        pass  # Optional - Quick pipelines


# =============================================================================
# COMPREHENSIVE PIPELINE METADATA
# =============================================================================

PIPELINE_METADATA = {
    # -------------------------------------------------------------------------
    # 1. HISAT2 + featureCounts
    # -------------------------------------------------------------------------
    "hisat2_featurecounts": {
        "module": 5,
        "stage": 2,
        "status": "implemented",
        "version": "2.2.0",
        "description": "HISAT2 alignment + featureCounts (low memory alternative)",
        "use_case": "Low memory systems (16GB), need BAM files",
        "memory_gb": 16,
        "produces_bam": True,
        "produces_gene_counts": True,
        "produces_isoform_counts": False,
        "supports_bootstraps": False,
        "validations": "15+",
        "typical_runtime": "30-60 min per sample",
        "requirements": ['hisat2>=2.2.0', 'featureCounts>=2.0.0'],
        "index_type": "HISAT2 index",
        "advantages": [
            "Lower memory than STAR (16GB vs 32GB)",
            "Produces BAM files for visualization",
            "Well-suited for gene-level analysis",
            "Fast alignment"
        ],
        "disadvantages": [
            "Slightly less accurate than STAR",
            "No isoform-level quantification",
            "No bootstrap support"
        ],
        "best_for": "Low memory systems (16GB), need BAM files"
    },
    
    # -------------------------------------------------------------------------
    # 2. Kallisto
    # -------------------------------------------------------------------------
    "kallisto": {
        "module": 5,
        "stage": 2,
        "status": "implemented",
        "version": "2.2.0",
        "description": "Kallisto pseudo-alignment (fastest, minimal memory)",
        "use_case": "Quick preliminary analysis, minimal resources",
        "memory_gb": 4,
        "produces_bam": False,
        "produces_gene_counts": True,
        "produces_isoform_counts": True,
        "supports_bootstraps": True,
        "validations": "12+",
        "typical_runtime": "5-10 min per sample",
        "requirements": ['kallisto>=0.48.0'],
        "index_type": "Kallisto transcriptome index",
        "advantages": [
            "Fastest pipeline (5-10 min/sample)",
            "Minimal memory (4GB)",
            "Gene + isoform quantification",
            "Bootstrap support for uncertainty"
        ],
        "disadvantages": [
            "No BAM files produced",
            "Slightly less accurate than alignment-based",
            "Requires transcriptome index"
        ],
        "best_for": "Quick preliminary analysis, minimal resources"
    },
    
    # -------------------------------------------------------------------------
    # 3. Salmon - RECOMMENDED
    # -------------------------------------------------------------------------
    "salmon": {
        "module": 5,
        "stage": 2,
        "status": "implemented",
        "version": "2.2.0",
        "description": "Salmon pseudo-alignment with bias correction",
        "use_case": "Standard differential expression (MOST POPULAR)",
        "memory_gb": 8,
        "produces_bam": False,
        "produces_gene_counts": True,
        "produces_isoform_counts": True,
        "supports_bootstraps": True,
        "validations": "15+",
        "typical_runtime": "10-20 min per sample",
        "requirements": ['salmon>=1.9.0'],
        "index_type": "Salmon decoy-aware index",
        "recommended": True,  # ⭐
        "advantages": [
            "Best balance speed/accuracy",
            "Sophisticated bias correction",
            "Gene + isoform quantification",
            "Bootstrap support",
            "Most popular for DE analysis"
        ],
        "disadvantages": [
            "No BAM files produced",
            "Requires transcriptome + decoy index"
        ],
        "best_for": "Standard differential expression (MOST POPULAR)"
    },
    
    # -------------------------------------------------------------------------
    # 4. STAR + featureCounts
    # -------------------------------------------------------------------------
    "star_featurecounts": {
        "module": 5,
        "stage": 2,
        "status": "implemented",
        "version": "2.2.0",
        "description": "STAR alignment + featureCounts (gold standard for genes)",
        "use_case": "Publication-quality gene-level analysis with BAM files",
        "memory_gb": 32,
        "produces_bam": True,
        "produces_gene_counts": True,
        "produces_isoform_counts": False,
        "supports_bootstraps": False,
        "validations": "17 (most comprehensive)",
        "typical_runtime": "40-70 min per sample",
        "requirements": ['STAR>=2.7.0', 'featureCounts>=2.0.0'],
        "index_type": "STAR genome index",
        "advantages": [
            "Gold standard for gene-level DE",
            "Produces BAM files",
            "Most comprehensive validations (17)",
            "Best for visualization (IGV)"
        ],
        "disadvantages": [
            "High memory (32GB)",
            "Slower than pseudo-alignment",
            "No isoform quantification",
            "No bootstrap support"
        ],
        "best_for": "Publication-quality gene-level analysis with BAM files"
    },
    
    # -------------------------------------------------------------------------
    # 5. STAR + RSEM
    # -------------------------------------------------------------------------
    "star_rsem": {
        "module": 5,
        "stage": 2,
        "status": "implemented",
        "version": "2.2.0",
        "description": "STAR + RSEM (gold standard for isoforms)",
        "use_case": "Isoform-level differential expression, alternative splicing",
        "memory_gb": 32,
        "produces_bam": True,
        "produces_gene_counts": True,
        "produces_isoform_counts": True,
        "supports_bootstraps": False,
        "supports_credibility_intervals": True,
        "validations": "15+",
        "typical_runtime": "60-120 min per sample",
        "requirements": ['STAR>=2.7.0', 'RSEM>=1.3.0'],
        "index_type": "STAR + RSEM transcriptome index",
        "advantages": [
            "Gold standard for isoform-level DE",
            "Gene + isoform quantification",
            "Produces BAM files",
            "Credibility intervals (optional)"
        ],
        "disadvantages": [
            "High memory (32GB)",
            "Slowest pipeline (60-120 min/sample)",
            "No bootstrap support"
        ],
        "best_for": "Isoform-level differential expression, alternative splicing"
    },
    
    # -------------------------------------------------------------------------
    # 6. STAR + Salmon - UNIQUE
    # -------------------------------------------------------------------------
    "star_salmon": {
        "module": 5,
        "stage": 2,
        "status": "implemented",
        "version": "2.2.0",
        "description": "STAR + Salmon hybrid (BAM + bootstraps)",
        "use_case": "Need BAM files AND bootstrap uncertainty (sleuth + IGV)",
        "memory_gb": 32,
        "produces_bam": True,
        "produces_gene_counts": True,
        "produces_isoform_counts": True,
        "supports_bootstraps": True,
        "validations": "15+",
        "typical_runtime": "50-90 min per sample",
        "requirements": ['STAR>=2.7.0', 'salmon>=1.9.0'],
        "index_type": "STAR genome index + Salmon transcriptome index",
        "unique": True,  # 🌟
        "advantages": [
            "UNIQUE: Only pipeline with BAM + bootstraps",
            "Best of both worlds",
            "Gene + isoform quantification",
            "Supports sleuth AND IGV"
        ],
        "disadvantages": [
            "High memory (32GB)",
            "Requires two indices",
            "Longer runtime"
        ],
        "best_for": "Need BAM files AND bootstrap uncertainty (sleuth + IGV)"
    },
    
    # -------------------------------------------------------------------------
    # 7. Quick Kallisto (Module 1)
    # -------------------------------------------------------------------------
    "quick_kallisto": {
        "module": 1,
        "stage": 1,
        "status": "planned",
        "version": "2.2.0",
        "description": "Fast Kallisto for QC profiling (no bootstraps)",
        "use_case": "Ultra-fast QC and preliminary profiling",
        "memory_gb": 4,
        "produces_bam": False,
        "produces_gene_counts": True,
        "produces_isoform_counts": True,
        "supports_bootstraps": False,
        "typical_runtime": "3-5 min per sample",
        "requirements": ['kallisto>=0.48.0'],
        "best_for": "Ultra-fast QC and preliminary profiling"
    },
    
    # -------------------------------------------------------------------------
    # 8. Quick Salmon (Module 1)
    # -------------------------------------------------------------------------
    "quick_salmon": {
        "module": 1,
        "stage": 1,
        "status": "planned",
        "version": "2.2.0",
        "description": "Fast Salmon for QC profiling (no bootstraps)",
        "use_case": "Fast QC with bias correction",
        "memory_gb": 8,
        "produces_bam": False,
        "produces_gene_counts": True,
        "produces_isoform_counts": True,
        "supports_bootstraps": False,
        "typical_runtime": "5-8 min per sample",
        "requirements": ['salmon>=1.9.0'],
        "best_for": "Fast QC with bias correction"
    }
}


# =============================================================================
# PIPELINE ACCESS FUNCTIONS
# =============================================================================

def get_pipeline(name: str) -> Optional[Type[BasePipeline]]:
    """
    Get pipeline class by name for instantiation and running.
    
    Parameters
    ----------
    name : str
        Pipeline name (e.g., 'salmon', 'star_featurecounts')
    
    Returns
    -------
    Type[BasePipeline] or None
        Pipeline class or None if not found
    
    Examples
    --------
    >>> SalmonPipeline = get_pipeline('salmon')
    >>> pipeline = SalmonPipeline(
    ...     output_dir='results/',
    ...     index='salmon_index',
    ...     bootstraps=100
    ... )
    >>> result = pipeline.run('samples.csv')
    """
    _load_pipelines()
    return _PIPELINE_REGISTRY.get(name)


def get_pipeline_params(name: str) -> Optional[Dict]:
    """
    Get CLI parameters for a pipeline.
    
    Parameters
    ----------
    name : str
        Pipeline name
    
    Returns
    -------
    Dict or None
        Parameter definitions for CLI integration
    
    Examples
    --------
    >>> params = get_pipeline_params('salmon')
    >>> print(params['bootstraps'])
    {'flag': '--bootstraps', 'default': 100, 'help': '...'}
    """
    _load_pipelines()
    return _PIPELINE_PARAMS.get(name)


def list_pipelines(module: int = None, status: str = None) -> List[str]:
    """
    List available pipeline names with optional filtering.
    
    Parameters
    ----------
    module : int, optional
        Filter by module number (1 or 5)
    status : str, optional
        Filter by status ('implemented' or 'planned')
    
    Returns
    -------
    List[str]
        Sorted list of pipeline names
    
    Examples
    --------
    >>> # List all pipelines
    >>> pipelines = list_pipelines()
    >>> print(pipelines)
    ['hisat2_featurecounts', 'kallisto', 'salmon', ...]
    
    >>> # List only implemented production pipelines
    >>> prod = list_pipelines(module=5, status='implemented')
    >>> print(f"Production pipelines: {len(prod)}")
    Production pipelines: 6
    """
    pipelines = []
    
    for name, info in PIPELINE_METADATA.items():
        # Apply filters
        if module is not None and info.get('module') != module:
            continue
        if status is not None and info.get('status') != status:
            continue
        pipelines.append(name)
    
    return sorted(pipelines)


def get_pipeline_info(name: str) -> Optional[Dict]:
    """
    Get comprehensive information about a pipeline.
    
    Parameters
    ----------
    name : str
        Pipeline name
    
    Returns
    -------
    Dict or None
        Pipeline metadata dictionary, or None if not found
    
    Examples
    --------
    >>> info = get_pipeline_info('salmon')
    >>> print(f"{info['description']}")
    Salmon pseudo-alignment with bias correction
    >>> print(f"Memory: {info['memory_gb']}GB")
    Memory: 8GB
    >>> print(f"Status: {info['status']}")
    Status: implemented
    """
    return PIPELINE_METADATA.get(name)


def get_implemented_pipelines() -> List[str]:
    """
    Get list of all implemented pipeline names.
    
    Returns
    -------
    List[str]
        List of implemented pipeline names
    
    Examples
    --------
    >>> implemented = get_implemented_pipelines()
    >>> print(f"✅ {len(implemented)} pipelines ready to use")
    ✅ 6 pipelines ready to use
    """
    return [name for name, info in PIPELINE_METADATA.items() 
            if info.get('status') == 'implemented']


# =============================================================================
# DISPLAY & COMPARISON FUNCTIONS
# =============================================================================

def print_pipeline_summary():
    """
    Print a formatted summary of all available pipelines.
    
    Examples
    --------
    >>> print_pipeline_summary()
    
    ══════════════════════════════════════════════════════════════════════════════════════════
      🦖 RAPTOR v2.2.0 - Pipeline Summary
    ══════════════════════════════════════════════════════════════════════════════════════════
    
      📦 Production Pipelines (Module 5): 6/6 Implemented ✅
      ──────────────────────────────────────────────────────────────────────────────────────
      ✅ hisat2_featurecounts  16GB   30-60 min       - HISAT2 alignment + featureCounts (low memory alternative)
      ✅ kallisto              4GB    5-10 min        - Kallisto pseudo-alignment (fastest, minimal memory)
      ✅ salmon ⭐             8GB    10-20 min       - Salmon pseudo-alignment with bias correction
      ...
    """
    print("\n" + "=" * 90)
    print("  🦖 RAPTOR v2.2.0 - Pipeline Summary")
    print("=" * 90)
    
    # Production pipelines
    production = [n for n, i in PIPELINE_METADATA.items() if i.get('module') == 5]
    implemented = [n for n in production if PIPELINE_METADATA[n].get('status') == 'implemented']
    
    print(f"\n  📦 Production Pipelines (Module 5): {len(implemented)}/6 Implemented ✅")
    print("  " + "-" * 86)
    
    for name in sorted(production):
        info = PIPELINE_METADATA[name]
        status = "✅" if info['status'] == 'implemented' else "⏳"
        mem = f"{info['memory_gb']}GB"
        time = info['typical_runtime'].split(' per')[0]
        
        # Add markers
        name_display = name
        if info.get('recommended'):
            name_display += " ⭐"
        if info.get('unique'):
            name_display += " 🌟"
        
        print(f"  {status} {name_display:<22} {mem:<6} {time:<15} - {info['description']}")
    
    # Quick pipelines
    quick = [n for n, i in PIPELINE_METADATA.items() if i.get('module') == 1]
    if quick:
        print(f"\n  ⚡ Quick Pipelines (Module 1): {len(quick)} Available")
        print("  " + "-" * 86)
        
        for name in sorted(quick):
            info = PIPELINE_METADATA[name]
            status = "⏳" if info['status'] == 'planned' else "✅"
            mem = f"{info['memory_gb']}GB"
            time = info['typical_runtime'].split(' per')[0]
            print(f"  {status} {name:<22} {mem:<6} {time:<15} - {info['description']}")
    
    print("\n" + "=" * 90)
    print(f"  Total Implemented: {len(get_implemented_pipelines())} pipelines")
    print("  ⭐ Recommended: salmon (best balance for most users)")
    print("  🌟 Unique: star_salmon (only pipeline with BAM + bootstraps)")
    print("=" * 90 + "\n")


def get_comparison_table() -> str:
    """
    Get a formatted comparison table of all production pipelines.
    
    Returns
    -------
    str
        Formatted comparison table
    
    Examples
    --------
    >>> print(get_comparison_table())
    
    ┌─────────────────────────────────────────────────────────────────────────────┐
    │                   🦖 RAPTOR v2.2.0 Pipeline Comparison                     │
    └─────────────────────────────────────────────────────────────────────────────┘
    
    Pipeline              BAM    Gene   Isoform   Bootstrap   Memory   Runtime
    ─────────────────────────────────────────────────────────────────────────────
    hisat2_featurecounts  ✅     ✅     ❌        ❌          16GB      30-60 min
    ...
    """
    lines = [
        "\n┌─────────────────────────────────────────────────────────────────────────────┐",
        "│                   🦖 RAPTOR v2.2.0 Pipeline Comparison                     │",
        "└─────────────────────────────────────────────────────────────────────────────┘",
        "",
        "Pipeline              BAM    Gene   Isoform   Bootstrap   Memory   Runtime",
        "─────────────────────────────────────────────────────────────────────────────"
    ]
    
    production = list_pipelines(module=5, status='implemented')
    for name in production:
        info = PIPELINE_METADATA[name]
        bam = "✅" if info.get('produces_bam') else "❌"
        gene = "✅" if info.get('produces_gene_counts') else "❌"
        iso = "✅" if info.get('produces_isoform_counts') else "❌"
        boot = "✅" if info.get('supports_bootstraps') else "❌"
        mem = f"{info['memory_gb']}GB"
        time = info['typical_runtime'].split(' per')[0]
        
        name_display = name[:20]
        if info.get('recommended'):
            name_display += " ⭐"
        if info.get('unique'):
            name_display += " 🌟"
        
        lines.append(f"{name_display:<22}{bam:<7}{gene:<7}{iso:<10}{boot:<12}{mem:<9}{time}")
    
    lines.append("─────────────────────────────────────────────────────────────────────────────")
    lines.append("⭐ Recommended   🌟 Unique capability")
    lines.append("")
    
    return '\n'.join(lines)


def get_pipeline_help(name: str) -> str:
    """
    Get detailed help text for a specific pipeline.
    
    Parameters
    ----------
    name : str
        Pipeline name
    
    Returns
    -------
    str
        Formatted help text with detailed information
    
    Examples
    --------
    >>> print(get_pipeline_help('salmon'))
    
    🦖 SALMON Pipeline ✅
    ══════════════════════════════════════════════════════════════════════════════════════════
    
    Version: 2.2.0
    Status: IMPLEMENTED
    Module: 5
    
    Salmon pseudo-alignment with bias correction
    
    🎯 Best Use Case:
       Standard differential expression (MOST POPULAR)
    ...
    """
    info = get_pipeline_info(name)
    if not info:
        available = ', '.join(get_implemented_pipelines())
        return f"❌ Pipeline '{name}' not found.\n\n✅ Available: {available}"
    
    status_emoji = "✅" if info['status'] == 'implemented' else "⏳"
    
    lines = [
        f"\n🦖 {name.upper().replace('_', ' + ')} Pipeline {status_emoji}",
        "=" * 90,
        f"\nVersion: {info.get('version', 'N/A')}",
        f"Status: {info['status'].upper()}",
        f"Module: {info['module']}",
        f"\n{info['description']}",
    ]
    
    if info.get('use_case'):
        lines.append(f"\n🎯 Best Use Case:\n   {info['use_case']}")
    
    lines.append(f"\n📊 Capabilities:")
    lines.append(f"   • Produces BAM files: {'Yes ✅' if info.get('produces_bam') else 'No'}")
    lines.append(f"   • Gene-level counts: {'Yes ✅' if info.get('produces_gene_counts') else 'No'}")
    lines.append(f"   • Isoform-level counts: {'Yes ✅' if info.get('produces_isoform_counts') else 'No'}")
    lines.append(f"   • Bootstrap support: {'Yes ✅' if info.get('supports_bootstraps') else 'No'}")
    if info.get('supports_credibility_intervals'):
        lines.append(f"   • Credibility intervals: Yes ✅ (optional)")
    
    lines.append(f"\n💻 Resources:")
    lines.append(f"   • Memory: ~{info['memory_gb']} GB")
    lines.append(f"   • Typical runtime: {info.get('typical_runtime', 'N/A')}")
    lines.append(f"   • Tools required: {', '.join(info.get('requirements', []))}")
    if info.get('index_type'):
        lines.append(f"   • Index type: {info['index_type']}")
    if info.get('validations'):
        lines.append(f"   • Validations: {info['validations']} checks")
    
    if info.get('advantages'):
        lines.append("\n✅ Advantages:")
        for adv in info['advantages']:
            lines.append(f"   • {adv}")
    
    if info.get('disadvantages'):
        lines.append("\n⚠️  Disadvantages:")
        for dis in info['disadvantages']:
            lines.append(f"   • {dis}")
    
    if info.get('best_for'):
        lines.append(f"\n💡 Best For: {info['best_for']}")
    
    lines.append("\n" + "=" * 90)
    
    return '\n'.join(lines)


def get_recommendation(requirements: Dict[str, bool]) -> str:
    """
    Get pipeline recommendation based on requirements.
    
    Parameters
    ----------
    requirements : dict
        Dictionary of requirements, e.g.:
        {
            'need_bam': True,
            'need_bootstraps': False,
            'low_memory': False,
            'need_isoforms': False
        }
    
    Returns
    -------
    str
        Recommended pipeline name with explanation
    
    Examples
    --------
    >>> rec = get_recommendation({'need_bam': True, 'need_bootstraps': True})
    >>> print(rec)
    'star_salmon (UNIQUE: only pipeline with BAM + bootstraps)'
    """
    need_bam = requirements.get('need_bam', False)
    need_bootstraps = requirements.get('need_bootstraps', False)
    low_memory = requirements.get('low_memory', False)
    need_isoforms = requirements.get('need_isoforms', False)
    speed_critical = requirements.get('speed_critical', False)
    
    # BAM + bootstraps = STAR + Salmon (unique!)
    if need_bam and need_bootstraps:
        return 'star_salmon (UNIQUE: only pipeline with BAM + bootstraps)'
    
    # BAM + isoforms = STAR + RSEM
    if need_bam and need_isoforms:
        return 'star_rsem (gold standard for isoforms with BAM files)'
    
    # BAM + genes only = STAR + featureCounts (or HISAT2 if low memory)
    if need_bam and not need_isoforms:
        if low_memory:
            return 'hisat2_featurecounts (low memory alternative with BAM)'
        return 'star_featurecounts (gold standard for gene-level with BAM)'
    
    # Speed critical + minimal memory = Kallisto
    if speed_critical:
        return 'kallisto (fastest, minimal memory, good for QC)'
    
    # Isoforms without BAM = Salmon (or Kallisto)
    if need_isoforms and not need_bam:
        if need_bootstraps:
            return 'salmon (best balance with bias correction and bootstraps)'
        return 'kallisto (faster than Salmon, good for preliminary analysis)'
    
    # Default: Salmon (most popular)
    return 'salmon (recommended for most users - best balance speed/accuracy)'


# =============================================================================
# MODULE EXPORTS
# =============================================================================

__all__ = [
    # Version
    '__version__',
    
    # Base classes
    'BasePipeline',
    'PipelineConfig',
    'PipelineResult',
    'SampleInfo',
    'ToolDependency',
    'DependencyHandler',
    
    # Utility functions
    'auto_detect_samples',
    'create_sample_sheet',
    
    # Pipeline access (for running)
    'get_pipeline',
    'get_pipeline_params',
    
    # Pipeline information (for display)
    'list_pipelines',
    'get_pipeline_info',
    'get_implemented_pipelines',
    
    # Display functions
    'print_pipeline_summary',
    'get_comparison_table',
    'get_pipeline_help',
    'get_recommendation',
    
    # Metadata
    'PIPELINE_METADATA'
]


# Auto-print summary when module is imported directly
if __name__ == '__main__':
    print_pipeline_summary()
    print(get_comparison_table())
