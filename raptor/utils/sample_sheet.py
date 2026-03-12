"""
RAPTOR v2.2.0 - Sample Sheet Utilities

Parse and validate RNA-seq sample sheets for pipelines.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
"""

import pandas as pd
from pathlib import Path
from typing import List, Optional, Union
from dataclasses import dataclass
import os

from .validation import (
    validate_file_path,
    validate_count_matrix
)
from .errors import ValidationError


@dataclass
class Sample:
    """
    Represents a single RNA-seq sample.
    
    Attributes
    ----------
    sample_id : str
        Unique sample identifier
    condition : str
        Experimental condition (e.g., 'Control', 'Treatment')
    batch : str
        Batch identifier for batch effect correction
    fastq_r1 : str
        Path to R1 FASTQ file
    fastq_r2 : str, optional
        Path to R2 FASTQ file (for paired-end)
    """
    sample_id: str
    condition: str
    batch: str
    fastq_r1: str
    fastq_r2: str = ''
    
    def __post_init__(self):
        """Validate sample after initialization."""
        if not self.sample_id:
            raise ValidationError("sample_id cannot be empty")
        
        if not self.fastq_r1:
            raise ValidationError(f"fastq_r1 missing for sample: {self.sample_id}")
    
    @property
    def is_paired(self) -> bool:
        """Check if sample is paired-end."""
        return bool(
            self.fastq_r2 and 
            self.fastq_r2 not in ['', 'NA', 'na', 'N/A', 'n/a', 'None', 'none']
        )
    
    @property
    def read_type(self) -> str:
        """Get read type ('paired' or 'single')."""
        return 'paired' if self.is_paired else 'single'
    
    def validate_files_exist(self) -> None:
        """
        Validate that FASTQ files exist.
        
        Raises
        ------
        ValidationError
            If any FASTQ file doesn't exist
        """
        if not os.path.exists(self.fastq_r1):
            raise ValidationError(
                f"FASTQ R1 not found for {self.sample_id}: {self.fastq_r1}"
            )
        
        if self.is_paired and not os.path.exists(self.fastq_r2):
            raise ValidationError(
                f"FASTQ R2 not found for {self.sample_id}: {self.fastq_r2}"
            )


class SampleSheet:
    """
    Parse and validate RNA-seq sample sheet.
    
    Sample sheet format (CSV):
    
    sample_id,condition,batch,fastq_r1,fastq_r2
    Sample1,Control,Batch1,/path/to/Sample1_R1.fastq.gz,/path/to/Sample1_R2.fastq.gz
    Sample2,Treatment,Batch1,/path/to/Sample2_R1.fastq.gz,/path/to/Sample2_R2.fastq.gz
    
    Parameters
    ----------
    filepath : str
        Path to sample sheet CSV
    validate_files : bool
        If True, validate that all FASTQ files exist
    
    Attributes
    ----------
    filepath : Path
        Path to sample sheet
    samples : List[Sample]
        List of Sample objects
    n_samples : int
        Number of samples
    
    Examples
    --------
    >>> from raptor.utils.sample_sheet import SampleSheet
    >>> sheet = SampleSheet('samples.csv')
    >>> print(f"Loaded {len(sheet.samples)} samples")
    >>> for sample in sheet.samples:
    ...     print(f"{sample.sample_id}: {sample.read_type}")
    """
    
    def __init__(self, filepath: Union[str, Path], validate_files: bool = False):
        """
        Initialize and load sample sheet.
        
        Parameters
        ----------
        filepath : Union[str, Path]
            Path to sample sheet CSV
        validate_files : bool
            If True, check that all FASTQ files exist
        """
        self.filepath = Path(filepath)
        self.samples: List[Sample] = []
        self.validate_files = validate_files
        
        self._load()
        self._validate()
    
    def _load(self) -> None:
        """Load sample sheet from CSV."""
        # Validate file exists
        if not self.filepath.exists():
            raise ValidationError(f"Sample sheet not found: {self.filepath}")
        
        # Load CSV
        try:
            df = pd.read_csv(self.filepath)
        except Exception as e:
            raise ValidationError(f"Failed to read sample sheet: {e}")
        
        # Check for empty
        if len(df) == 0:
            raise ValidationError("Sample sheet is empty")
        
        # Validate required columns
        required_columns = ['sample_id', 'fastq_r1']
        missing_columns = [col for col in required_columns if col not in df.columns]
        
        if missing_columns:
            raise ValidationError(
                f"Sample sheet missing required columns: {', '.join(missing_columns)}\n"
                f"Required: {', '.join(required_columns)}\n"
                f"Found: {', '.join(df.columns)}"
            )
        
        # Optional columns (provide defaults if missing)
        if 'condition' not in df.columns:
            df['condition'] = ''
        if 'batch' not in df.columns:
            df['batch'] = ''
        if 'fastq_r2' not in df.columns:
            df['fastq_r2'] = ''
        
        # Convert NaN to empty string
        df = df.fillna('')
        
        # Create Sample objects
        for _, row in df.iterrows():
            try:
                sample = Sample(
                    sample_id=str(row['sample_id']),
                    condition=str(row.get('condition', '')),
                    batch=str(row.get('batch', '')),
                    fastq_r1=str(row['fastq_r1']),
                    fastq_r2=str(row.get('fastq_r2', ''))
                )
                self.samples.append(sample)
            except Exception as e:
                raise ValidationError(f"Error parsing row {_ + 1}: {e}")
    
    def _validate(self) -> None:
        """Validate sample sheet."""
        if not self.samples:
            raise ValidationError("No samples loaded from sample sheet")
        
        # Check for duplicate sample IDs
        sample_ids = [s.sample_id for s in self.samples]
        duplicates = [sid for sid in set(sample_ids) if sample_ids.count(sid) > 1]
        
        if duplicates:
            raise ValidationError(
                f"Duplicate sample IDs found: {', '.join(duplicates)}"
            )
        
        # Validate FASTQ files exist (if requested)
        if self.validate_files:
            for sample in self.samples:
                sample.validate_files_exist()
    
    @property
    def n_samples(self) -> int:
        """Get number of samples."""
        return len(self.samples)
    
    @property
    def read_types(self) -> List[str]:
        """Get list of read types for all samples."""
        return [s.read_type for s in self.samples]
    
    @property
    def is_mixed_reads(self) -> bool:
        """Check if sample sheet contains both paired and single-end reads."""
        types = set(self.read_types)
        return len(types) > 1
    
    @property
    def conditions(self) -> List[str]:
        """Get list of unique conditions."""
        return sorted(set(s.condition for s in self.samples if s.condition))
    
    @property
    def batches(self) -> List[str]:
        """Get list of unique batches."""
        return sorted(set(s.batch for s in self.samples if s.batch))
    
    def get_samples_by_condition(self, condition: str) -> List[Sample]:
        """
        Get all samples for a specific condition.
        
        Parameters
        ----------
        condition : str
            Condition name
        
        Returns
        -------
        List[Sample]
            Samples matching condition
        """
        return [s for s in self.samples if s.condition == condition]
    
    def get_samples_by_batch(self, batch: str) -> List[Sample]:
        """
        Get all samples for a specific batch.
        
        Parameters
        ----------
        batch : str
            Batch identifier
        
        Returns
        -------
        List[Sample]
            Samples matching batch
        """
        return [s for s in self.samples if s.batch == batch]
    
    def to_dataframe(self) -> pd.DataFrame:
        """
        Convert sample sheet to DataFrame.
        
        Returns
        -------
        pd.DataFrame
            Sample sheet as DataFrame
        """
        data = []
        for sample in self.samples:
            data.append({
                'sample_id': sample.sample_id,
                'condition': sample.condition,
                'batch': sample.batch,
                'fastq_r1': sample.fastq_r1,
                'fastq_r2': sample.fastq_r2,
                'read_type': sample.read_type
            })
        return pd.DataFrame(data)
    
    def save(self, output_path: Union[str, Path]) -> None:
        """
        Save sample sheet to CSV.
        
        Parameters
        ----------
        output_path : Union[str, Path]
            Output file path
        """
        df = self.to_dataframe()
        df.to_csv(output_path, index=False)
    
    def __len__(self) -> int:
        """Get number of samples."""
        return len(self.samples)
    
    def __iter__(self):
        """Iterate over samples."""
        return iter(self.samples)
    
    def __getitem__(self, index: int) -> Sample:
        """Get sample by index."""
        return self.samples[index]
    
    def __repr__(self) -> str:
        """String representation."""
        return (
            f"SampleSheet(n_samples={self.n_samples}, "
            f"conditions={len(self.conditions)}, "
            f"batches={len(self.batches)})"
        )


def load_sample_sheet(
    filepath: Union[str, Path],
    validate_files: bool = False
) -> SampleSheet:
    """
    Convenience function to load sample sheet.
    
    Parameters
    ----------
    filepath : Union[str, Path]
        Path to sample sheet CSV
    validate_files : bool
        If True, validate that all FASTQ files exist
    
    Returns
    -------
    SampleSheet
        Loaded sample sheet
    
    Examples
    --------
    >>> from raptor.utils.sample_sheet import load_sample_sheet
    >>> sheet = load_sample_sheet('samples.csv', validate_files=True)
    >>> print(f"Loaded {len(sheet)} samples")
    """
    return SampleSheet(filepath, validate_files=validate_files)


__all__ = [
    'Sample',
    'SampleSheet',
    'load_sample_sheet',
]
