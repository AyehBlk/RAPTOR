#!/usr/bin/env python3

"""
Sample Sheet Handler for RAPTOR Pipelines

Handles parsing and validation of sample sheets for both paired-end
and single-end FASTQ files. Used by all RAPTOR pipelines.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0

Example Usage:
    >>> from raptor import SampleSheet
    >>> sheet = SampleSheet('samples.csv')
    >>> sheet.validate()
    >>> for sample in sheet:
    ...     print(sample.sample_id, sample.fastq_r1)

Sample Sheet Format (Paired-end):
    sample_id,condition,batch,fastq_r1,fastq_r2
    Sample1,Control,Batch1,Sample1_R1.fastq.gz,Sample1_R2.fastq.gz
    Sample2,Treatment,Batch1,Sample2_R1.fastq.gz,Sample2_R2.fastq.gz

Sample Sheet Format (Single-end):
    sample_id,condition,batch,fastq
    Sample1,Control,Batch1,Sample1.fastq.gz
    Sample2,Treatment,Batch1,Sample2.fastq.gz
"""

import os
import re
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass

import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class Sample:
    """
    Represents a single sample with FASTQ file(s).
    
    Parameters
    ----------
    sample_id : str
        Unique identifier for the sample
    fastq_r1 : str
        Path to R1 FASTQ file (or single-end file)
    fastq_r2 : str, optional
        Path to R2 FASTQ file (None for single-end)
    condition : str, optional
        Experimental condition (e.g., Control, Treatment)
    batch : str, optional
        Batch identifier for batch effect correction
    
    Attributes
    ----------
    is_paired : bool
        True if sample has paired-end reads
    
    Examples
    --------
    >>> sample = Sample(
    ...     sample_id='Sample1',
    ...     fastq_r1='Sample1_R1.fastq.gz',
    ...     fastq_r2='Sample1_R2.fastq.gz',
    ...     condition='Control'
    ... )
    >>> sample.is_paired
    True
    """
    sample_id: str
    fastq_r1: str
    fastq_r2: Optional[str] = None
    condition: Optional[str] = None
    batch: Optional[str] = None
    
    @property
    def is_paired(self) -> bool:
        """Check if sample is paired-end."""
        return self.fastq_r2 is not None
    
    def validate(self) -> Tuple[bool, List[str]]:
        """
        Validate sample files exist.
        
        Returns
        -------
        Tuple[bool, List[str]]
            (is_valid, list of error messages)
        """
        errors = []
        
        # Check R1 exists
        if not os.path.exists(self.fastq_r1):
            errors.append(f"R1 file not found: {self.fastq_r1}")
        
        # Check R2 exists (if paired)
        if self.fastq_r2 and not os.path.exists(self.fastq_r2):
            errors.append(f"R2 file not found: {self.fastq_r2}")
        
        # Check file extensions
        valid_extensions = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']
        r1_valid = any(self.fastq_r1.endswith(ext) for ext in valid_extensions)
        if not r1_valid:
            errors.append(f"Invalid R1 file extension: {self.fastq_r1}")
        
        if self.fastq_r2:
            r2_valid = any(self.fastq_r2.endswith(ext) for ext in valid_extensions)
            if not r2_valid:
                errors.append(f"Invalid R2 file extension: {self.fastq_r2}")
        
        return len(errors) == 0, errors


class SampleSheet:
    """
    Parse and validate sample sheets for RNA-seq analysis.
    
    Supports two formats:
    
    1. Paired-end:
       sample_id,fastq_r1,fastq_r2[,condition,batch,...]
       
    2. Single-end:
       sample_id,fastq[,condition,batch,...]
    
    Parameters
    ----------
    filepath : str or Path
        Path to sample sheet CSV/TSV file
    
    Attributes
    ----------
    samples : List[Sample]
        List of parsed samples
    is_paired : bool
        Whether samples are paired-end
    metadata_columns : List[str]
        Additional metadata columns found
    
    Examples
    --------
    >>> sheet = SampleSheet('samples.csv')
    >>> sheet.validate()
    (True, [])
    >>> for sample in sheet.samples:
    ...     print(sample.sample_id, sample.fastq_r1)
    Sample1 /path/to/Sample1_R1.fastq.gz
    Sample2 /path/to/Sample2_R1.fastq.gz
    
    >>> # Access by index
    >>> sheet[0].sample_id
    'Sample1'
    
    >>> # Number of samples
    >>> len(sheet)
    6
    """
    
    # Required columns for each type
    REQUIRED_PAIRED = ['sample_id', 'fastq_r1', 'fastq_r2']
    REQUIRED_SINGLE = ['sample_id', 'fastq']
    
    # Alternative column names (flexible parsing)
    COLUMN_ALIASES = {
        'sample_id': ['sample', 'sample_name', 'name', 'id', 'sampleid', 'sample_id'],
        'fastq_r1': ['fastq_r1', 'fastq_1', 'r1', 'read1', 'fq1', 'fastq1'],
        'fastq_r2': ['fastq_r2', 'fastq_2', 'r2', 'read2', 'fq2', 'fastq2'],
        'fastq': ['fastq', 'fq', 'reads', 'read', 'file'],
        'condition': ['condition', 'group', 'treatment', 'class'],
        'batch': ['batch', 'batch_id', 'run', 'lane']
    }
    
    def __init__(self, filepath: str):
        """Initialize sample sheet from file."""
        self.filepath = Path(filepath)
        self.samples: List[Sample] = []
        self.is_paired: bool = False
        self.metadata_columns: List[str] = []
        self._df: Optional[pd.DataFrame] = None
        
        self._load()
    
    def _load(self):
        """Load and parse sample sheet file."""
        if not self.filepath.exists():
            raise FileNotFoundError(f"Sample sheet not found: {self.filepath}")
        
        # Detect separator
        suffix = self.filepath.suffix.lower()
        if suffix in ['.tsv', '.txt']:
            sep = '\t'
        else:
            sep = ','
        
        # Load file
        self._df = pd.read_csv(self.filepath, sep=sep, dtype=str)
        self._df.columns = self._df.columns.str.strip().str.lower()
        
        # Normalize column names
        self._normalize_columns()
        
        # Detect paired vs single
        self._detect_read_type()
        
        # Parse samples
        self._parse_samples()
        
        logger.info(f"Loaded {len(self.samples)} samples ({self._read_type})")
    
    def _normalize_columns(self):
        """Normalize column names using aliases."""
        column_map = {}
        
        for standard_name, aliases in self.COLUMN_ALIASES.items():
            for alias in aliases:
                if alias in self._df.columns:
                    column_map[alias] = standard_name
                    break
        
        self._df = self._df.rename(columns=column_map)
    
    def _detect_read_type(self):
        """Detect if samples are paired-end or single-end."""
        has_r1 = 'fastq_r1' in self._df.columns
        has_r2 = 'fastq_r2' in self._df.columns
        has_single = 'fastq' in self._df.columns
        
        if has_r1 and has_r2:
            self.is_paired = True
            self._read_type = 'paired-end'
        elif has_single:
            self.is_paired = False
            self._read_type = 'single-end'
        elif has_r1 and not has_r2:
            # Only R1 provided - treat as single-end
            self.is_paired = False
            self._read_type = 'single-end'
            self._df['fastq'] = self._df['fastq_r1']
        else:
            raise ValueError(
                "Cannot determine read type. Sample sheet must have either:\n"
                "  - 'fastq_r1' and 'fastq_r2' columns (paired-end)\n"
                "  - 'fastq' column (single-end)"
            )
        
        # Identify metadata columns
        required = set(['sample_id', 'fastq_r1', 'fastq_r2', 'fastq'])
        self.metadata_columns = [c for c in self._df.columns if c not in required]
    
    def _parse_samples(self):
        """Parse samples from dataframe."""
        self.samples = []
        
        for _, row in self._df.iterrows():
            sample_id = str(row['sample_id']).strip()
            
            if self.is_paired:
                fastq_r1 = str(row['fastq_r1']).strip()
                fastq_r2 = str(row['fastq_r2']).strip()
            else:
                fastq_r1 = str(row['fastq']).strip()
                fastq_r2 = None
            
            # Get optional metadata
            condition = row.get('condition', None)
            if pd.notna(condition):
                condition = str(condition).strip()
            else:
                condition = None
            
            batch = row.get('batch', None)
            if pd.notna(batch):
                batch = str(batch).strip()
            else:
                batch = None
            
            sample = Sample(
                sample_id=sample_id,
                fastq_r1=fastq_r1,
                fastq_r2=fastq_r2,
                condition=condition,
                batch=batch
            )
            
            self.samples.append(sample)
    
    def validate(self, check_files: bool = True) -> Tuple[bool, List[str]]:
        """
        Validate sample sheet.
        
        Parameters
        ----------
        check_files : bool
            Whether to check if FASTQ files exist
        
        Returns
        -------
        Tuple[bool, List[str]]
            (is_valid, list of error messages)
        """
        errors = []
        
        # Check for empty sample sheet
        if len(self.samples) == 0:
            errors.append("Sample sheet is empty")
            return False, errors
        
        # Check for duplicate sample IDs
        sample_ids = [s.sample_id for s in self.samples]
        duplicates = set([x for x in sample_ids if sample_ids.count(x) > 1])
        if duplicates:
            errors.append(f"Duplicate sample IDs: {duplicates}")
        
        # Validate each sample
        if check_files:
            for sample in self.samples:
                is_valid, sample_errors = sample.validate()
                if not is_valid:
                    errors.extend(sample_errors)
        
        return len(errors) == 0, errors
    
    def get_sample_ids(self) -> List[str]:
        """Get list of sample IDs."""
        return [s.sample_id for s in self.samples]
    
    def get_sample(self, sample_id: str) -> Optional[Sample]:
        """
        Get sample by ID.
        
        Parameters
        ----------
        sample_id : str
            Sample ID to find
        
        Returns
        -------
        Sample or None
            Sample object if found, None otherwise
        """
        for s in self.samples:
            if s.sample_id == sample_id:
                return s
        return None
    
    def get_conditions(self) -> Optional[Dict[str, str]]:
        """Get sample-to-condition mapping if available."""
        if not any(s.condition for s in self.samples):
            return None
        return {s.sample_id: s.condition for s in self.samples}
    
    def get_batches(self) -> Optional[Dict[str, str]]:
        """Get sample-to-batch mapping if available."""
        if not any(s.batch for s in self.samples):
            return None
        return {s.sample_id: s.batch for s in self.samples}
    
    def to_dataframe(self) -> pd.DataFrame:
        """
        Convert to pandas DataFrame.
        
        Returns
        -------
        pd.DataFrame
            DataFrame with all sample information
        """
        return self._df.copy()
    
    def to_metadata_df(self) -> pd.DataFrame:
        """
        Convert to metadata DataFrame for downstream analysis.
        
        Returns
        -------
        pd.DataFrame
            DataFrame with sample_id as index
        """
        data = []
        for s in self.samples:
            row = {
                'sample_id': s.sample_id,
                'read_type': 'paired' if s.is_paired else 'single'
            }
            if s.condition:
                row['condition'] = s.condition
            if s.batch:
                row['batch'] = s.batch
            data.append(row)
        
        return pd.DataFrame(data).set_index('sample_id')
    
    def summary(self) -> str:
        """Generate human-readable summary."""
        lines = [
            "=" * 60,
            "SAMPLE SHEET SUMMARY",
            "=" * 60,
            f"File: {self.filepath}",
            f"Read type: {self._read_type.upper()}",
            f"Samples: {len(self.samples)}",
            ""
        ]
        
        # Conditions summary
        conditions = self.get_conditions()
        if conditions:
            cond_counts = {}
            for cond in conditions.values():
                cond_counts[cond] = cond_counts.get(cond, 0) + 1
            lines.append("Conditions:")
            for cond, count in cond_counts.items():
                lines.append(f"  - {cond}: {count} samples")
            lines.append("")
        
        # Sample table
        lines.append(f"{'Sample ID':<20} {'Condition':<15} {'R1':<25}")
        lines.append("-" * 60)
        
        for s in self.samples[:10]:  # Show first 10
            cond = s.condition or "-"
            r1_name = Path(s.fastq_r1).name[:23] if s.fastq_r1 else "-"
            lines.append(f"{s.sample_id:<20} {cond:<15} {r1_name:<25}")
        
        if len(self.samples) > 10:
            lines.append(f"... and {len(self.samples) - 10} more samples")
        
        lines.append("=" * 60)
        
        return "\n".join(lines)
    
    def __len__(self) -> int:
        return len(self.samples)
    
    def __iter__(self):
        return iter(self.samples)
    
    def __getitem__(self, index):
        return self.samples[index]
    
    def __repr__(self):
        return f"SampleSheet({self.filepath}, n={len(self.samples)}, {self._read_type})"


def auto_detect_samples(input_dir: str, 
                        pattern: str = 'auto') -> Tuple[List[Sample], str]:
    """
    Auto-detect FASTQ files and create samples.
    
    Scans a directory for FASTQ files and automatically pairs R1/R2 files.
    
    Parameters
    ----------
    input_dir : str
        Directory containing FASTQ files
    pattern : str
        Naming pattern or 'auto' for automatic detection
    
    Returns
    -------
    Tuple[List[Sample], str]
        (list of samples, read type 'paired' or 'single')
    
    Examples
    --------
    >>> samples, read_type = auto_detect_samples('fastq/')
    >>> print(f"Found {len(samples)} {read_type} samples")
    Found 6 paired samples
    
    >>> for s in samples:
    ...     print(s.sample_id)
    Sample1
    Sample2
    """
    input_path = Path(input_dir)
    
    if not input_path.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    
    # Find all FASTQ files
    fastq_patterns = ['*.fastq.gz', '*.fq.gz', '*.fastq', '*.fq']
    fastq_files = []
    for pat in fastq_patterns:
        fastq_files.extend(input_path.glob(pat))
    
    if not fastq_files:
        raise ValueError(f"No FASTQ files found in: {input_dir}")
    
    fastq_files = sorted([str(f) for f in fastq_files])
    
    # R1 indicators and their R2 replacements
    r1_to_r2 = [
        ('_R1_001', '_R2_001'),  # Illumina: sample_R1_001.fastq.gz
        ('_R1.', '_R2.'),        # Standard: sample_R1.fastq.gz
        ('_R1_', '_R2_'),        # With suffix: sample_R1_trimmed.fastq.gz
        ('.R1.', '.R2.'),        # Dot notation: sample.R1.fastq.gz
        ('_1.', '_2.'),          # Simple: sample_1.fastq.gz
    ]
    
    samples = []
    used_files = set()
    
    # Try to find pairs
    for fq in fastq_files:
        if fq in used_files:
            continue
        
        filename = os.path.basename(fq)
        found_pair = False
        
        for r1_marker, r2_marker in r1_to_r2:
            if r1_marker in filename:
                # This looks like R1, find R2
                r2_filename = filename.replace(r1_marker, r2_marker)
                r2_path = os.path.join(os.path.dirname(fq), r2_filename)
                
                if r2_path in fastq_files:
                    # Found pair!
                    # Extract sample name (remove R1 marker and extension)
                    sample_name = filename.split(r1_marker)[0]
                    # Clean up sample name
                    sample_name = re.sub(r'[_\-\.]$', '', sample_name)
                    
                    samples.append(Sample(
                        sample_id=sample_name,
                        fastq_r1=fq,
                        fastq_r2=r2_path
                    ))
                    used_files.add(fq)
                    used_files.add(r2_path)
                    found_pair = True
                    break
        
        # If no pair found, treat as single-end
        if not found_pair and fq not in used_files:
            # Extract sample name
            sample_name = os.path.basename(fq)
            for ext in ['.fastq.gz', '.fq.gz', '.fastq', '.fq']:
                sample_name = sample_name.replace(ext, '')
            # Remove R1/R2 suffixes if present
            sample_name = re.sub(r'[_\.]?(R?[12]|r?[12])$', '', sample_name)
            
            samples.append(Sample(
                sample_id=sample_name,
                fastq_r1=fq,
                fastq_r2=None
            ))
            used_files.add(fq)
    
    # Determine read type
    paired_count = sum(1 for s in samples if s.is_paired)
    single_count = len(samples) - paired_count
    
    if paired_count > 0 and single_count > 0:
        logger.warning(f"Mixed read types detected: {paired_count} paired, {single_count} single")
    
    read_type = 'paired' if paired_count >= single_count else 'single'
    
    logger.info(f"Auto-detected {len(samples)} {read_type} samples")
    
    return samples, read_type


def create_sample_sheet_template(output_path: str, 
                                  paired: bool = True,
                                  n_samples: int = 6) -> str:
    """
    Create a sample sheet template file.
    
    Parameters
    ----------
    output_path : str
        Path to save template
    paired : bool
        Whether to create paired-end template
    n_samples : int
        Number of example samples
    
    Returns
    -------
    str
        Path to created template
    
    Examples
    --------
    >>> create_sample_sheet_template('samples.csv', paired=True)
    'samples.csv'
    
    >>> # Then edit the file with your actual paths
    """
    if paired:
        header = "sample_id,condition,batch,fastq_r1,fastq_r2"
        rows = []
        for i in range(1, n_samples + 1):
            cond = "Control" if i <= n_samples // 2 else "Treatment"
            batch = f"Batch{(i-1) % 2 + 1}"
            rows.append(
                f"Sample{i},{cond},{batch},"
                f"/path/to/Sample{i}_R1.fastq.gz,"
                f"/path/to/Sample{i}_R2.fastq.gz"
            )
    else:
        header = "sample_id,condition,batch,fastq"
        rows = []
        for i in range(1, n_samples + 1):
            cond = "Control" if i <= n_samples // 2 else "Treatment"
            batch = f"Batch{(i-1) % 2 + 1}"
            rows.append(
                f"Sample{i},{cond},{batch},"
                f"/path/to/Sample{i}.fastq.gz"
            )
    
    content = header + "\n" + "\n".join(rows) + "\n"
    
    with open(output_path, 'w') as f:
        f.write(content)
    
    logger.info(f"Created sample sheet template: {output_path}")
    
    return output_path


if __name__ == '__main__':
    # Example usage
    print("RAPTOR Sample Sheet Handler")
    print("=" * 40)
    print()
    print("Create a sample sheet template:")
    print("  from raptor import create_sample_sheet_template")
    print("  create_sample_sheet_template('samples.csv', paired=True)")
    print()
    print("Load and validate a sample sheet:")
    print("  from raptor import SampleSheet")
    print("  sheet = SampleSheet('samples.csv')")
    print("  is_valid, errors = sheet.validate()")
    print("  print(sheet.summary())")
    print()
    print("Auto-detect samples from directory:")
    print("  from raptor import auto_detect_samples")
    print("  samples, read_type = auto_detect_samples('fastq/')")
