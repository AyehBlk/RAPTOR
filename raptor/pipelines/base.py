#!/usr/bin/env python3
"""
RAPTOR v2.2.0 - Base Pipeline Classes

Provides base classes and dependency handling for all production pipelines.

Dependency Handling Strategy (Hybrid - Option F):
=================================================
1. Default: Check PATH, error with clear instructions if missing
2. Flag --use-docker: Use containerized version
3. Flag --modules: Load HPC modules before execution

Author: Ayeh Bolouki
Email: ayehbolouki1988@gmail.com
Version: 2.2.0
"""

import os
import sys
import logging
import subprocess
import shutil
import json
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Any, Union, Tuple
from datetime import datetime
import pandas as pd

logger = logging.getLogger(__name__)


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class ToolDependency:
    """
    Defines a tool dependency with installation options.
    
    Attributes
    ----------
    name : str
        Human-readable name (e.g., "STAR Aligner")
    command : str
        Command to check/run (e.g., "STAR")
    min_version : str
        Minimum required version (e.g., "2.7.0")
    docker_image : str
        Docker image containing the tool
    conda_package : str
        Conda package name for installation
    hpc_module : str
        HPC module name (e.g., "STAR/2.7.10b")
    install_url : str
        URL for manual installation instructions
    install_conda : str
        Conda install command
    version_flag : str
        Flag to get version (default: --version)
    version_pattern : str
        Regex pattern to extract version from output
    """
    name: str
    command: str
    min_version: str = ""
    docker_image: str = ""
    conda_package: str = ""
    hpc_module: str = ""
    install_url: str = ""
    install_conda: str = ""
    version_flag: str = "--version"
    version_pattern: str = r"(\d+\.\d+\.?\d*)"


@dataclass
class SampleInfo:
    """
    Information about a single sample.
    
    Attributes
    ----------
    sample_id : str
        Unique sample identifier
    fastq_1 : Path
        Path to R1 FASTQ file (or single-end file)
    fastq_2 : Path or None
        Path to R2 FASTQ file (None for single-end)
    condition : str
        Experimental condition/group
    batch : str
        Batch identifier (optional)
    extra : Dict
        Additional metadata
    """
    sample_id: str
    fastq_1: Path
    fastq_2: Optional[Path] = None
    condition: str = ""
    batch: str = ""
    extra: Dict[str, Any] = field(default_factory=dict)
    
    @property
    def is_paired(self) -> bool:
        """Check if sample is paired-end."""
        return self.fastq_2 is not None
    
    def validate(self) -> bool:
        """Validate that FASTQ files exist."""
        if not self.fastq_1.exists():
            raise FileNotFoundError(f"FASTQ R1 not found: {self.fastq_1}")
        if self.fastq_2 and not self.fastq_2.exists():
            raise FileNotFoundError(f"FASTQ R2 not found: {self.fastq_2}")
        return True
    
    @classmethod
    def from_dict(cls, data: Dict) -> 'SampleInfo':
        """Create SampleInfo from dictionary."""
        return cls(
            sample_id=data['sample_id'],
            fastq_1=Path(data.get('fastq_r1') or data.get('fastq_1') or data.get('fastq')),
            fastq_2=Path(data['fastq_r2']) if data.get('fastq_r2') or data.get('fastq_2') else None,
            condition=data.get('condition', ''),
            batch=data.get('batch', ''),
            extra={k: v for k, v in data.items() 
                   if k not in ['sample_id', 'fastq_r1', 'fastq_r2', 'fastq_1', 'fastq_2', 
                               'fastq', 'condition', 'batch']}
        )


@dataclass
class PipelineConfig:
    """
    Pipeline configuration.
    
    Attributes
    ----------
    output_dir : Path
        Main output directory
    threads : int
        Number of threads to use
    memory_gb : int
        Memory limit in GB
    use_docker : bool
        Whether to use Docker execution
    docker_image : str
        Docker image to use
    modules : List[str]
        HPC modules to load
    keep_bam : bool
        Whether to keep BAM files (for pipelines that produce them)
    extra_args : Dict
        Additional tool-specific arguments
    """
    output_dir: Path
    threads: int = 8
    memory_gb: int = 32
    use_docker: bool = False
    docker_image: str = ""
    modules: List[str] = field(default_factory=list)
    keep_bam: bool = True
    extra_args: Dict[str, Any] = field(default_factory=dict)


@dataclass
class PipelineResult:
    """
    Result from pipeline execution.
    
    Attributes
    ----------
    success : bool
        Whether pipeline completed successfully
    gene_counts_file : Path or None
        Path to gene-level count matrix
    tx_counts_file : Path or None
        Path to transcript-level count matrix (if applicable)
    tpm_file : Path or None
        Path to TPM matrix
    bam_dir : Path or None
        Directory containing BAM files (if applicable)
    logs_dir : Path
        Directory containing log files
    sample_results : List[Dict]
        Per-sample results
    errors : List[str]
        Any errors encountered
    warnings : List[str]
        Any warnings
    start_time : datetime
        Pipeline start time
    end_time : datetime
        Pipeline end time
    metrics : Dict
        Execution metrics (timing, resource usage)
    """
    success: bool
    gene_counts_file: Optional[Path] = None
    tx_counts_file: Optional[Path] = None
    tpm_file: Optional[Path] = None
    bam_dir: Optional[Path] = None
    logs_dir: Optional[Path] = None
    sample_results: List[Dict] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None
    metrics: Dict[str, Any] = field(default_factory=dict)
    
    @property
    def duration_minutes(self) -> float:
        """Get pipeline duration in minutes."""
        if self.start_time and self.end_time:
            return (self.end_time - self.start_time).total_seconds() / 60
        return 0.0
    
    def to_dict(self) -> Dict:
        """Convert result to dictionary."""
        return {
            'success': self.success,
            'gene_counts_file': str(self.gene_counts_file) if self.gene_counts_file else None,
            'tx_counts_file': str(self.tx_counts_file) if self.tx_counts_file else None,
            'tpm_file': str(self.tpm_file) if self.tpm_file else None,
            'bam_dir': str(self.bam_dir) if self.bam_dir else None,
            'logs_dir': str(self.logs_dir) if self.logs_dir else None,
            'n_samples': len(self.sample_results),
            'errors': self.errors,
            'warnings': self.warnings,
            'duration_minutes': self.duration_minutes,
            'metrics': self.metrics
        }
    
    def summary(self) -> str:
        """Generate human-readable summary."""
        lines = []
        status = "✅ SUCCESS" if self.success else "❌ FAILED"
        lines.append(f"\n{status}")
        lines.append("=" * 50)
        
        if self.gene_counts_file:
            lines.append(f"📊 Gene counts: {self.gene_counts_file}")
        if self.tx_counts_file:
            lines.append(f"📊 Tx counts: {self.tx_counts_file}")
        if self.tpm_file:
            lines.append(f"📊 TPM: {self.tpm_file}")
        if self.bam_dir:
            lines.append(f"🗂️  BAM files: {self.bam_dir}")
        
        lines.append(f"📋 Samples processed: {len(self.sample_results)}")
        lines.append(f"⏱️  Duration: {self.duration_minutes:.1f} min")
        
        if self.errors:
            lines.append(f"\n❌ Errors:")
            for err in self.errors:
                lines.append(f"   • {err}")
        
        if self.warnings:
            lines.append(f"\n⚠️  Warnings:")
            for warn in self.warnings:
                lines.append(f"   • {warn}")
        
        return '\n'.join(lines)


# =============================================================================
# DEPENDENCY HANDLER (Hybrid Approach)
# =============================================================================

class DependencyHandler:
    """
    Handles tool dependencies with multiple resolution strategies.
    
    Strategy (Hybrid - Option F):
    =============================
    1. Default: Check if tool exists in PATH, error with instructions if not
    2. --use-docker: Use Docker container with tool pre-installed
    3. --modules: Load HPC environment modules before execution
    
    Examples
    --------
    >>> handler = DependencyHandler()
    >>> # Check if STAR is available
    >>> handler.check_dependency(star_dep, use_docker=False, modules=[])
    >>> # Get execution prefix (docker run / module load / empty)
    >>> prefix = handler.get_execution_prefix(use_docker=True, docker_image='star:2.7.10b')
    """
    
    def __init__(self):
        self._docker_available = None
        self._singularity_available = None
    
    def check_dependency(
        self,
        dep: ToolDependency,
        use_docker: bool = False,
        docker_image: str = "",
        modules: List[str] = None
    ) -> Tuple[bool, str]:
        """
        Check if a dependency is available.
        
        Parameters
        ----------
        dep : ToolDependency
            Dependency to check
        use_docker : bool
            If True, assume Docker will provide the tool
        docker_image : str
            Docker image to use
        modules : List[str]
            HPC modules that might provide the tool
        
        Returns
        -------
        Tuple[bool, str]
            (available, message)
        """
        modules = modules or []
        
        # If Docker mode, assume available via container
        if use_docker:
            image = docker_image or dep.docker_image
            if not image:
                return False, f"No Docker image specified for {dep.name}"
            
            if not self._check_docker_available():
                return False, "Docker not available. Install Docker or use --modules instead."
            
            return True, f"{dep.name} will be provided by Docker image: {image}"
        
        # If modules specified, check if they're loadable
        if modules:
            # In HPC environment, we trust the modules
            return True, f"{dep.name} will be loaded via modules: {', '.join(modules)}"
        
        # Default: Check PATH
        return self._check_in_path(dep)
    
    def _check_in_path(self, dep: ToolDependency) -> Tuple[bool, str]:
        """Check if tool is available in PATH."""
        tool_path = shutil.which(dep.command)
        
        if not tool_path:
            # Generate helpful error message
            msg = self._generate_install_instructions(dep)
            return False, msg
        
        # Check version if specified
        if dep.min_version:
            version = self._get_tool_version(dep)
            if version and self._compare_versions(version, dep.min_version) < 0:
                return False, (
                    f"{dep.name} version {version} found, but {dep.min_version}+ required.\n"
                    f"Update: {dep.install_conda or dep.install_url}"
                )
        
        return True, f"{dep.name} found at: {tool_path}"
    
    def _generate_install_instructions(self, dep: ToolDependency) -> str:
        """Generate helpful installation instructions."""
        lines = [
            f"\n❌ {dep.name} not found in PATH",
            "",
            "📋 Installation options:",
        ]
        
        if dep.install_conda:
            lines.append(f"   • Conda: {dep.install_conda}")
        
        if dep.docker_image:
            lines.append(f"   • Docker: raptor pipeline --use-docker --docker-image {dep.docker_image}")
        
        if dep.hpc_module:
            lines.append(f"   • HPC: raptor pipeline --modules \"{dep.hpc_module}\"")
        
        if dep.install_url:
            lines.append(f"   • Manual: {dep.install_url}")
        
        return '\n'.join(lines)
    
    def _get_tool_version(self, dep: ToolDependency) -> Optional[str]:
        """Get tool version."""
        import re
        try:
            result = subprocess.run(
                [dep.command, dep.version_flag],
                capture_output=True,
                text=True,
                timeout=30
            )
            output = result.stdout + result.stderr
            match = re.search(dep.version_pattern, output)
            return match.group(1) if match else None
        except Exception:
            return None
    
    def _compare_versions(self, v1: str, v2: str) -> int:
        """Compare two version strings. Returns -1 if v1 < v2, 0 if equal, 1 if v1 > v2."""
        def normalize(v):
            return [int(x) for x in re.sub(r'(\.0+)*$', '', v).split('.')]
        import re
        try:
            n1, n2 = normalize(v1), normalize(v2)
            return (n1 > n2) - (n1 < n2)
        except:
            return 0
    
    def _check_docker_available(self) -> bool:
        """Check if Docker is available."""
        if self._docker_available is None:
            self._docker_available = shutil.which('docker') is not None
        return self._docker_available
    
    def _check_singularity_available(self) -> bool:
        """Check if Singularity is available."""
        if self._singularity_available is None:
            self._singularity_available = (
                shutil.which('singularity') is not None or
                shutil.which('apptainer') is not None
            )
        return self._singularity_available
    
    def get_execution_wrapper(
        self,
        use_docker: bool = False,
        docker_image: str = "",
        modules: List[str] = None,
        volumes: Dict[str, str] = None
    ) -> 'ExecutionWrapper':
        """
        Get an execution wrapper for running commands.
        
        Parameters
        ----------
        use_docker : bool
            Use Docker execution
        docker_image : str
            Docker image to use
        modules : List[str]
            HPC modules to load
        volumes : Dict[str, str]
            Volume mappings for Docker (host:container)
        
        Returns
        -------
        ExecutionWrapper
            Wrapper object for command execution
        """
        modules = modules or []
        volumes = volumes or {}
        
        if use_docker:
            return DockerExecutionWrapper(docker_image, volumes)
        elif modules:
            return ModuleExecutionWrapper(modules)
        else:
            return LocalExecutionWrapper()


# =============================================================================
# EXECUTION WRAPPERS
# =============================================================================

class ExecutionWrapper(ABC):
    """Abstract base class for execution wrappers."""
    
    @abstractmethod
    def run(
        self,
        cmd: List[str],
        cwd: Optional[Path] = None,
        env: Optional[Dict[str, str]] = None,
        log_file: Optional[Path] = None
    ) -> Tuple[int, str, str]:
        """
        Run a command.
        
        Returns
        -------
        Tuple[int, str, str]
            (returncode, stdout, stderr)
        """
        pass
    
    @abstractmethod
    def get_prefix(self) -> List[str]:
        """Get command prefix (e.g., 'docker run ...')."""
        pass


class LocalExecutionWrapper(ExecutionWrapper):
    """Execute commands directly on the local system."""
    
    def run(
        self,
        cmd: List[str],
        cwd: Optional[Path] = None,
        env: Optional[Dict[str, str]] = None,
        log_file: Optional[Path] = None
    ) -> Tuple[int, str, str]:
        """Run command locally."""
        try:
            full_env = os.environ.copy()
            if env:
                full_env.update(env)
            
            result = subprocess.run(
                cmd,
                cwd=cwd,
                env=full_env,
                capture_output=True,
                text=True
            )
            
            if log_file:
                with open(log_file, 'w') as f:
                    f.write(f"Command: {' '.join(cmd)}\n")
                    f.write(f"Return code: {result.returncode}\n")
                    f.write("\n=== STDOUT ===\n")
                    f.write(result.stdout)
                    f.write("\n=== STDERR ===\n")
                    f.write(result.stderr)
            
            return result.returncode, result.stdout, result.stderr
            
        except Exception as e:
            error_msg = str(e)
            if log_file:
                with open(log_file, 'w') as f:
                    f.write(f"Error: {error_msg}\n")
            return 1, "", error_msg
    
    def get_prefix(self) -> List[str]:
        """No prefix for local execution."""
        return []


class DockerExecutionWrapper(ExecutionWrapper):
    """Execute commands in a Docker container."""
    
    def __init__(self, image: str, volumes: Dict[str, str] = None):
        self.image = image
        self.volumes = volumes or {}
    
    def run(
        self,
        cmd: List[str],
        cwd: Optional[Path] = None,
        env: Optional[Dict[str, str]] = None,
        log_file: Optional[Path] = None
    ) -> Tuple[int, str, str]:
        """Run command in Docker container."""
        docker_cmd = ['docker', 'run', '--rm']
        
        # Add volume mounts
        for host_path, container_path in self.volumes.items():
            docker_cmd.extend(['-v', f'{host_path}:{container_path}'])
        
        # Auto-mount working directory
        if cwd:
            docker_cmd.extend(['-v', f'{cwd}:/work', '-w', '/work'])
        
        # Add environment variables
        if env:
            for key, value in env.items():
                docker_cmd.extend(['-e', f'{key}={value}'])
        
        # Add image and command
        docker_cmd.append(self.image)
        docker_cmd.extend(cmd)
        
        # Execute
        local_wrapper = LocalExecutionWrapper()
        return local_wrapper.run(docker_cmd, cwd=None, env=None, log_file=log_file)
    
    def get_prefix(self) -> List[str]:
        """Get Docker prefix."""
        prefix = ['docker', 'run', '--rm']
        for host_path, container_path in self.volumes.items():
            prefix.extend(['-v', f'{host_path}:{container_path}'])
        prefix.append(self.image)
        return prefix


class ModuleExecutionWrapper(ExecutionWrapper):
    """Execute commands with HPC modules loaded."""
    
    def __init__(self, modules: List[str]):
        self.modules = modules
    
    def run(
        self,
        cmd: List[str],
        cwd: Optional[Path] = None,
        env: Optional[Dict[str, str]] = None,
        log_file: Optional[Path] = None
    ) -> Tuple[int, str, str]:
        """Run command with modules loaded."""
        # Build shell command with module loads
        module_loads = ' && '.join([f'module load {m}' for m in self.modules])
        full_cmd = f"{module_loads} && {' '.join(cmd)}"
        
        shell_cmd = ['bash', '-c', full_cmd]
        
        local_wrapper = LocalExecutionWrapper()
        return local_wrapper.run(shell_cmd, cwd=cwd, env=env, log_file=log_file)
    
    def get_prefix(self) -> List[str]:
        """No prefix, modules are loaded in shell."""
        return []


# =============================================================================
# BASE PIPELINE CLASS
# =============================================================================

class BasePipeline(ABC):
    """
    Abstract base class for all production pipelines.
    
    Implements common functionality:
    - Dependency checking with Hybrid strategy
    - Sample sheet parsing
    - Output directory setup
    - Logging
    - Results aggregation
    
    Subclasses must implement:
    - get_dependencies(): Return list of tool dependencies
    - run_sample(): Process a single sample
    - combine_counts(): Aggregate counts from all samples
    
    Parameters
    ----------
    output_dir : str or Path
        Output directory for results
    threads : int
        Number of threads (default: 8)
    memory_gb : int
        Memory limit in GB (default: 32)
    use_docker : bool
        Use Docker for execution (default: False)
    docker_image : str
        Docker image to use (default: pipeline's default)
    modules : List[str]
        HPC modules to load (default: None)
    keep_bam : bool
        Keep BAM files if produced (default: True)
    extra_args : Dict
        Additional tool-specific arguments
    
    Examples
    --------
    >>> # Local execution (tools in PATH)
    >>> pipeline = SalmonPipeline(output_dir='results/', threads=8)
    >>> 
    >>> # Docker execution
    >>> pipeline = StarFeatureCountsPipeline(
    ...     output_dir='results/',
    ...     use_docker=True,
    ...     docker_image='quay.io/biocontainers/star:2.7.10b'
    ... )
    >>> 
    >>> # HPC with modules
    >>> pipeline = StarRsemPipeline(
    ...     output_dir='results/',
    ...     modules=['STAR/2.7.10b', 'RSEM/1.3.3']
    ... )
    """
    
    # Override in subclasses
    PIPELINE_NAME: str = "base"
    PIPELINE_VERSION: str = "2.2.0"
    DESCRIPTION: str = "Base pipeline"
    DOCKER_IMAGE: str = ""
    
    def __init__(
        self,
        output_dir: Union[str, Path],
        threads: int = 8,
        memory_gb: int = 32,
        use_docker: bool = False,
        docker_image: Optional[str] = None,
        modules: Optional[List[str]] = None,
        keep_bam: bool = True,
        extra_args: Optional[Dict[str, str]] = None
    ):
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.memory_gb = memory_gb
        self.use_docker = use_docker
        self.docker_image = docker_image or self.DOCKER_IMAGE
        self.modules = modules or []
        self.keep_bam = keep_bam
        self.extra_args = extra_args or {}
        
        # Setup directories
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.logs_dir = self.output_dir / "logs"
        self.logs_dir.mkdir(exist_ok=True)
        
        # Initialize handlers
        self.dep_handler = DependencyHandler()
        self._execution_wrapper = None
    
    @abstractmethod
    def get_dependencies(self) -> List[ToolDependency]:
        """
        Return list of required dependencies.
        
        Returns
        -------
        List[ToolDependency]
            List of tool dependencies
        """
        pass
    
    @abstractmethod
    def run_sample(
        self,
        sample: SampleInfo,
        index_path: Path,
        **kwargs
    ) -> Dict[str, Any]:
        """
        Process a single sample.
        
        Parameters
        ----------
        sample : SampleInfo
            Sample information
        index_path : Path
            Path to reference index
        **kwargs
            Additional arguments
        
        Returns
        -------
        Dict
            Result dictionary with at least:
            - 'success': bool
            - 'sample_id': str
            - 'counts_file': str (if successful)
            - 'error': str (if failed)
        """
        pass
    
    @abstractmethod
    def combine_counts(
        self,
        sample_results: List[Dict],
        output_file: Path
    ) -> pd.DataFrame:
        """
        Combine per-sample counts into a matrix.
        
        Parameters
        ----------
        sample_results : List[Dict]
            Results from run_sample for each sample
        output_file : Path
            Output file path
        
        Returns
        -------
        pd.DataFrame
            Gene count matrix (genes × samples)
        """
        pass
    
    def check_dependencies(self) -> Tuple[bool, List[str]]:
        """
        Check all pipeline dependencies.
        
        Returns
        -------
        Tuple[bool, List[str]]
            (all_ok, messages)
        """
        deps = self.get_dependencies()
        all_ok = True
        messages = []
        
        for dep in deps:
            ok, msg = self.dep_handler.check_dependency(
                dep,
                use_docker=self.use_docker,
                docker_image=self.docker_image,
                modules=self.modules
            )
            messages.append(msg)
            if not ok:
                all_ok = False
        
        return all_ok, messages
    
    def get_execution_wrapper(self) -> ExecutionWrapper:
        """Get the execution wrapper based on configuration."""
        if self._execution_wrapper is None:
            # Setup volumes for Docker
            volumes = {}
            if self.use_docker:
                # Mount output directory
                volumes[str(self.output_dir.absolute())] = '/output'
            
            self._execution_wrapper = self.dep_handler.get_execution_wrapper(
                use_docker=self.use_docker,
                docker_image=self.docker_image,
                modules=self.modules,
                volumes=volumes
            )
        return self._execution_wrapper
    
    def run_command(
        self,
        cmd: List[str],
        log_file: Optional[Path] = None,
        cwd: Optional[Path] = None,
        env: Optional[Dict[str, str]] = None
    ) -> Tuple[int, str, str]:
        """
        Run a command using the configured execution wrapper.
        
        Parameters
        ----------
        cmd : List[str]
            Command and arguments
        log_file : Path, optional
            File to write logs
        cwd : Path, optional
            Working directory
        env : Dict[str, str], optional
            Environment variables
        
        Returns
        -------
        Tuple[int, str, str]
            (returncode, stdout, stderr)
        """
        wrapper = self.get_execution_wrapper()
        return wrapper.run(cmd, cwd=cwd, env=env, log_file=log_file)
    
    def load_samples(
        self,
        sample_sheet: Union[str, Path, pd.DataFrame]
    ) -> List[SampleInfo]:
        """
        Load samples from sample sheet.
        
        Parameters
        ----------
        sample_sheet : str, Path, or DataFrame
            Sample sheet file or DataFrame
        
        Returns
        -------
        List[SampleInfo]
            List of sample information objects
        """
        if isinstance(sample_sheet, pd.DataFrame):
            df = sample_sheet
        else:
            df = pd.read_csv(sample_sheet)
        
        samples = []
        for _, row in df.iterrows():
            sample = SampleInfo.from_dict(row.to_dict())
            samples.append(sample)
        
        return samples
    
    def run(
        self,
        sample_sheet: Union[str, Path, pd.DataFrame],
        index_path: Union[str, Path],
        **kwargs
    ) -> PipelineResult:
        """
        Run the complete pipeline.
        
        Parameters
        ----------
        sample_sheet : str, Path, or DataFrame
            Sample sheet with sample information
        index_path : str or Path
            Path to reference index
        **kwargs
            Additional pipeline-specific arguments
        
        Returns
        -------
        PipelineResult
            Pipeline execution result
        """
        start_time = datetime.now()
        index_path = Path(index_path)
        
        logger.info(f"🦖 Starting {self.PIPELINE_NAME} pipeline v{self.PIPELINE_VERSION}")
        logger.info(f"   Output: {self.output_dir}")
        logger.info(f"   Threads: {self.threads}")
        
        # Check dependencies
        logger.info("📋 Checking dependencies...")
        deps_ok, dep_messages = self.check_dependencies()
        for msg in dep_messages:
            logger.info(f"   {msg}")
        
        if not deps_ok:
            return PipelineResult(
                success=False,
                errors=["Dependency check failed"] + dep_messages,
                start_time=start_time,
                end_time=datetime.now()
            )
        
        # Load samples
        logger.info("📂 Loading samples...")
        try:
            samples = self.load_samples(sample_sheet)
            logger.info(f"   Found {len(samples)} samples")
        except Exception as e:
            return PipelineResult(
                success=False,
                errors=[f"Failed to load samples: {e}"],
                start_time=start_time,
                end_time=datetime.now()
            )
        
        # Process samples
        logger.info("🚀 Processing samples...")
        sample_results = []
        errors = []
        warnings = []
        
        for i, sample in enumerate(samples, 1):
            logger.info(f"   [{i}/{len(samples)}] {sample.sample_id}")
            try:
                sample.validate()
                result = self.run_sample(sample, index_path, **kwargs)
                sample_results.append(result)
                
                if not result.get('success', False):
                    errors.append(f"{sample.sample_id}: {result.get('error', 'Unknown error')}")
                    
            except Exception as e:
                errors.append(f"{sample.sample_id}: {str(e)}")
                sample_results.append({
                    'success': False,
                    'sample_id': sample.sample_id,
                    'error': str(e)
                })
        
        # Check for failures
        successful_results = [r for r in sample_results if r.get('success', False)]
        if not successful_results:
            return PipelineResult(
                success=False,
                errors=errors,
                sample_results=sample_results,
                logs_dir=self.logs_dir,
                start_time=start_time,
                end_time=datetime.now()
            )
        
        # Combine counts
        logger.info("📊 Combining counts...")
        gene_counts_file = self.output_dir / "gene_counts.csv"
        try:
            counts_df = self.combine_counts(successful_results, gene_counts_file)
            logger.info(f"   Created: {gene_counts_file}")
            logger.info(f"   Shape: {counts_df.shape[0]} genes × {counts_df.shape[1]} samples")
        except Exception as e:
            errors.append(f"Failed to combine counts: {e}")
            return PipelineResult(
                success=False,
                errors=errors,
                sample_results=sample_results,
                logs_dir=self.logs_dir,
                start_time=start_time,
                end_time=datetime.now()
            )
        
        # Generate TPM if method exists
        tpm_file = None
        if hasattr(self, 'generate_tpm'):
            try:
                tpm_file = self.output_dir / "tpm.csv"
                self.generate_tpm(successful_results, tpm_file)
                logger.info(f"   Created: {tpm_file}")
            except Exception as e:
                warnings.append(f"TPM generation failed: {e}")
        
        # Determine BAM directory
        bam_dir = None
        if self.keep_bam:
            bam_path = self.output_dir / "bam"
            if bam_path.exists() and any(bam_path.glob('*.bam')):
                bam_dir = bam_path
        
        end_time = datetime.now()
        
        # Save run info
        run_info = {
            'pipeline': self.PIPELINE_NAME,
            'version': self.PIPELINE_VERSION,
            'start_time': start_time.isoformat(),
            'end_time': end_time.isoformat(),
            'duration_minutes': (end_time - start_time).total_seconds() / 60,
            'threads': self.threads,
            'n_samples': len(samples),
            'n_successful': len(successful_results),
            'use_docker': self.use_docker,
            'docker_image': self.docker_image if self.use_docker else None,
            'modules': self.modules if self.modules else None,
            'extra_args': self.extra_args
        }
        with open(self.output_dir / 'pipeline_info.json', 'w') as f:
            json.dump(run_info, f, indent=2)
        
        logger.info(f"\n✅ Pipeline complete in {run_info['duration_minutes']:.1f} minutes")
        
        return PipelineResult(
            success=True,
            gene_counts_file=gene_counts_file,
            tpm_file=tpm_file,
            bam_dir=bam_dir,
            logs_dir=self.logs_dir,
            sample_results=sample_results,
            errors=errors,
            warnings=warnings,
            start_time=start_time,
            end_time=end_time,
            metrics=run_info
        )


# =============================================================================
# SAMPLE SHEET UTILITIES
# =============================================================================

def auto_detect_samples(
    fastq_dir: Union[str, Path],
    pattern: str = None
) -> Tuple[List[SampleInfo], str]:
    """
    Auto-detect samples from a FASTQ directory.
    
    Parameters
    ----------
    fastq_dir : str or Path
        Directory containing FASTQ files
    pattern : str, optional
        Glob pattern for FASTQ files
    
    Returns
    -------
    Tuple[List[SampleInfo], str]
        (samples, read_type) where read_type is 'paired' or 'single'
    """
    fastq_dir = Path(fastq_dir)
    
    # Common FASTQ patterns
    patterns = pattern.split(',') if pattern else [
        '*_R1*.fastq.gz', '*_R1*.fq.gz',
        '*_1.fastq.gz', '*_1.fq.gz',
        '*.fastq.gz', '*.fq.gz'
    ]
    
    # Find R1 files
    r1_files = []
    for pat in patterns:
        r1_files.extend(fastq_dir.glob(pat))
    r1_files = list(set(r1_files))
    
    if not r1_files:
        raise FileNotFoundError(f"No FASTQ files found in {fastq_dir}")
    
    samples = []
    read_type = 'unknown'
    
    for r1_path in sorted(r1_files):
        # Try to find R2
        r1_name = r1_path.name
        r2_path = None
        
        # Common R1 -> R2 patterns
        r2_patterns = [
            (r'_R1', '_R2'),
            (r'_R1_', '_R2_'),
            (r'_1\.', '_2.'),
            (r'\.R1\.', '.R2.')
        ]
        
        import re
        for r1_pat, r2_pat in r2_patterns:
            r2_name = re.sub(r1_pat, r2_pat, r1_name)
            if r2_name != r1_name:
                potential_r2 = r1_path.parent / r2_name
                if potential_r2.exists():
                    r2_path = potential_r2
                    break
        
        # Extract sample ID
        sample_id = r1_name
        for suffix in ['.fastq.gz', '.fq.gz', '.fastq', '.fq']:
            sample_id = sample_id.replace(suffix, '')
        for pat in ['_R1', '_R1_001', '_1', '.R1']:
            sample_id = sample_id.replace(pat, '')
        
        samples.append(SampleInfo(
            sample_id=sample_id,
            fastq_1=r1_path,
            fastq_2=r2_path
        ))
        
        if r2_path:
            read_type = 'paired'
        elif read_type == 'unknown':
            read_type = 'single'
    
    return samples, read_type


def create_sample_sheet(
    samples: List[SampleInfo],
    output_file: Union[str, Path]
) -> pd.DataFrame:
    """
    Create a sample sheet CSV from sample info objects.
    
    Parameters
    ----------
    samples : List[SampleInfo]
        List of sample information
    output_file : str or Path
        Output CSV file path
    
    Returns
    -------
    pd.DataFrame
        Sample sheet DataFrame
    """
    data = []
    for s in samples:
        row = {
            'sample_id': s.sample_id,
            'condition': s.condition,
            'batch': s.batch,
            'fastq_r1': str(s.fastq_1),
            'fastq_r2': str(s.fastq_2) if s.fastq_2 else ''
        }
        row.update(s.extra)
        data.append(row)
    
    df = pd.DataFrame(data)
    df.to_csv(output_file, index=False)
    return df


# =============================================================================
# EXPORTS
# =============================================================================

__all__ = [
    'ToolDependency',
    'SampleInfo',
    'PipelineConfig',
    'PipelineResult',
    'DependencyHandler',
    'ExecutionWrapper',
    'LocalExecutionWrapper',
    'DockerExecutionWrapper',
    'ModuleExecutionWrapper',
    'BasePipeline',
    'auto_detect_samples',
    'create_sample_sheet'
]
