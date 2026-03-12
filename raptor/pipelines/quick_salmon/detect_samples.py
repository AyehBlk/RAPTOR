#!/usr/bin/env python3

"""
Auto-detect FASTQ samples from directory.

Creates sample sheet CSV for quantification pipelines.
Shared utility for both Salmon and Kallisto pipelines.

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
"""

import argparse
import sys
import os
import glob
import re
from pathlib import Path
from typing import List, Tuple, Optional
from dataclasses import dataclass
import pandas as pd


@dataclass
class DetectedSample:
    """Represents a detected FASTQ sample."""
    sample_id: str
    fastq_r1: str
    fastq_r2: Optional[str] = None
    
    @property
    def is_paired(self) -> bool:
        return self.fastq_r2 is not None


class FastqAutoDetector:
    """
    Auto-detect FASTQ files and pair R1/R2.
    
    Supports common naming conventions:
    - Illumina: Sample_S1_L001_R1_001.fastq.gz
    - Generic: Sample_R1.fastq.gz, Sample_1.fastq.gz
    """
    
    FASTQ_EXTENSIONS = ['.fastq.gz', '.fq.gz', '.fastq', '.fq']
    
    R1_PATTERNS = [
        r'_R1_001\.', r'_R1\.', r'\.R1\.', r'_1\.',
        r'_read1\.', r'_r1\.', r'\.1\.'
    ]
    
    R2_PATTERNS = [
        r'_R2_001\.', r'_R2\.', r'\.R2\.', r'_2\.',
        r'_read2\.', r'_r2\.', r'\.2\.'
    ]
    
    def __init__(self, input_dir: str, recursive: bool = False):
        """
        Initialize detector.
        
        Parameters
        ----------
        input_dir : str
            Directory containing FASTQ files
        recursive : bool
            Search subdirectories
        """
        self.input_dir = Path(input_dir)
        self.recursive = recursive
    
    def find_fastq_files(self) -> List[str]:
        """Find all FASTQ files in directory."""
        files = []
        
        for ext in self.FASTQ_EXTENSIONS:
            if self.recursive:
                pattern = str(self.input_dir / '**' / f'*{ext}')
                files.extend(glob.glob(pattern, recursive=True))
            else:
                pattern = str(self.input_dir / f'*{ext}')
                files.extend(glob.glob(pattern))
        
        return sorted(set(files))
    
    def classify_files(self, files: List[str]) -> Tuple[List[str], List[str], List[str]]:
        """
        Classify files as R1, R2, or unknown.
        
        Returns
        -------
        Tuple[List[str], List[str], List[str]]
            (r1_files, r2_files, unknown_files)
        """
        r1_files = []
        r2_files = []
        unknown_files = []
        
        for f in files:
            filename = os.path.basename(f)
            
            is_r1 = any(re.search(p, filename, re.IGNORECASE) for p in self.R1_PATTERNS)
            is_r2 = any(re.search(p, filename, re.IGNORECASE) for p in self.R2_PATTERNS)
            
            if is_r1 and not is_r2:
                r1_files.append(f)
            elif is_r2 and not is_r1:
                r2_files.append(f)
            else:
                unknown_files.append(f)
        
        return r1_files, r2_files, unknown_files
    
    def extract_sample_name(self, filepath: str) -> str:
        """Extract sample name from FASTQ filename."""
        filename = os.path.basename(filepath)
        
        # Remove extension
        for ext in self.FASTQ_EXTENSIONS:
            if filename.endswith(ext):
                filename = filename[:-len(ext)]
                break
        
        # Remove read identifiers
        patterns_to_remove = [
            r'_R1_001$', r'_R2_001$',
            r'_R1$', r'_R2$',
            r'\.R1$', r'\.R2$',
            r'_1$', r'_2$',
            r'_read1$', r'_read2$',
            r'_r1$', r'_r2$'
        ]
        
        for pattern in patterns_to_remove:
            filename = re.sub(pattern, '', filename, flags=re.IGNORECASE)
        
        # Clean up trailing underscores/dots
        filename = re.sub(r'[_\.]+$', '', filename)
        
        return filename
    
    def pair_samples(self, r1_files: List[str], r2_files: List[str]) -> List[DetectedSample]:
        """Pair R1 and R2 files by sample name."""
        r1_by_sample = {self.extract_sample_name(f): f for f in r1_files}
        r2_by_sample = {self.extract_sample_name(f): f for f in r2_files}
        
        samples = []
        all_sample_names = set(r1_by_sample.keys()) | set(r2_by_sample.keys())
        
        for sample_name in sorted(all_sample_names):
            r1 = r1_by_sample.get(sample_name)
            r2 = r2_by_sample.get(sample_name)
            
            if r1 and r2:
                samples.append(DetectedSample(sample_name, r1, r2))
            elif r1:
                print(f"  ⚠️  {sample_name}: R1 found but no matching R2")
                samples.append(DetectedSample(sample_name, r1, None))
            elif r2:
                print(f"  ⚠️  {sample_name}: R2 found but no matching R1 (skipping)")
        
        return samples
    
    def detect(self) -> Tuple[List[DetectedSample], str]:
        """
        Detect samples from directory.
        
        Returns
        -------
        Tuple[List[DetectedSample], str]
            (samples, read_type)
        """
        # Find all FASTQ files
        files = self.find_fastq_files()
        
        if not files:
            return [], 'unknown'
        
        # Classify files
        r1_files, r2_files, unknown_files = self.classify_files(files)
        
        print(f"📂 Found {len(files)} FASTQ files:")
        print(f"   R1 files: {len(r1_files)}")
        print(f"   R2 files: {len(r2_files)}")
        if unknown_files:
            print(f"   Unclassified: {len(unknown_files)}")
        
        # Determine read type
        if r2_files:
            read_type = 'paired'
            samples = self.pair_samples(r1_files, r2_files)
        else:
            read_type = 'single'
            # Use R1 files or unknown files for single-end
            files_to_use = r1_files if r1_files else unknown_files
            samples = [
                DetectedSample(self.extract_sample_name(f), f, None)
                for f in files_to_use
            ]
        
        print(f"   Read type: {read_type}")
        print(f"   Samples: {len(samples)}")
        
        return samples, read_type


def create_sample_sheet(
    samples: List[DetectedSample],
    output_path: str,
    read_type: str = 'paired'
) -> pd.DataFrame:
    """
    Create sample sheet CSV from detected samples.
    
    Parameters
    ----------
    samples : List[DetectedSample]
        Detected samples
    output_path : str
        Output CSV path
    read_type : str
        'paired' or 'single'
    
    Returns
    -------
    pd.DataFrame
        Sample sheet dataframe
    """
    if read_type == 'paired':
        data = [{
            'sample_id': s.sample_id,
            'condition': '',
            'batch': '',
            'fastq_r1': s.fastq_r1,
            'fastq_r2': s.fastq_r2 or ''
        } for s in samples]
    else:
        data = [{
            'sample_id': s.sample_id,
            'condition': '',
            'batch': '',
            'fastq_r1': s.fastq_r1,
            'fastq_r2': ''
        } for s in samples]
    
    df = pd.DataFrame(data)
    df.to_csv(output_path, index=False)
    
    return df


def auto_detect_samples(input_dir: str) -> Tuple[List[DetectedSample], str]:
    """
    Convenience function for auto-detecting samples.
    
    This function is used by the CLI and other modules.
    
    Parameters
    ----------
    input_dir : str
        Directory containing FASTQ files
    
    Returns
    -------
    Tuple[List[DetectedSample], str]
        (samples, read_type)
    """
    detector = FastqAutoDetector(input_dir)
    return detector.detect()


def main():
    parser = argparse.ArgumentParser(
        description='Auto-detect FASTQ files and create sample sheet (RAPTOR v2.2.0)'
    )
    
    parser.add_argument(
        '--input', '-i',
        required=True,
        help='Directory containing FASTQ files'
    )
    
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output CSV file path'
    )
    
    parser.add_argument(
        '--read-type', '-r',
        choices=['auto', 'paired', 'single'],
        default='auto',
        help='Read type (default: auto-detect)'
    )
    
    parser.add_argument(
        '--recursive',
        action='store_true',
        help='Search subdirectories'
    )
    
    args = parser.parse_args()
    
    # Validate input directory
    if not os.path.isdir(args.input):
        print(f"❌ Error: Directory not found: {args.input}")
        sys.exit(1)
    
    # Detect samples
    detector = FastqAutoDetector(args.input, recursive=args.recursive)
    samples, detected_read_type = detector.detect()
    
    if not samples:
        print("❌ No FASTQ files found!")
        sys.exit(1)
    
    # Use specified read type or detected
    read_type = detected_read_type if args.read_type == 'auto' else args.read_type
    
    # Create sample sheet
    create_sample_sheet(samples, args.output, read_type)
    
    print(f"\n✅ Sample sheet created: {args.output}")
    print(f"   {len(samples)} samples ({read_type}-end)")
    print(f"\n📝 Next steps:")
    print(f"   1. Edit {args.output} to add condition and batch columns")
    print(f"   2. Run: raptor quick-count -m salmon -s {args.output} -i <index>")


if __name__ == '__main__':
    main()
