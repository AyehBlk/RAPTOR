"""
Real-Time Resource Monitoring

Tracks CPU, memory, disk I/O, and other system resources during pipeline execution
to optimize performance and identify bottlenecks.

Author: Ayeh Bolouki
Email: ayeh.bolouki@unamur.be
"""

import os
import time
import threading
import json
from pathlib import Path
from typing import Dict, Optional, List
import logging

import psutil
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

logger = logging.getLogger(__name__)


class ResourceMonitor:
    """
    Monitor system resources during analysis.
    
    Tracks CPU usage, memory consumption, disk I/O, and network activity
    in real-time during pipeline execution.
    
    Parameters
    ----------
    output_dir : str
        Directory for monitoring outputs
    interval : float
        Sampling interval in seconds (default: 1.0)
    max_memory_mb : int, optional
        Maximum memory threshold in MB for alerts
    max_cpu_percent : int, optional
        Maximum CPU usage threshold for alerts
    
    Examples
    --------
    >>> monitor = ResourceMonitor(output_dir="monitoring")
    >>> monitor.start_monitoring()
    >>> # Run your analysis here
    >>> stats = monitor.stop_monitoring()
    >>> print(f"Peak memory: {stats['memory']['peak_mb']:.1f} MB")
    >>> monitor.generate_plots()
    """
    
    def __init__(
        self,
        output_dir: str = "resource_monitoring",
        interval: float = 1.0,
        max_memory_mb: Optional[int] = None,
        max_cpu_percent: Optional[int] = None
    ):
        """Initialize resource monitor."""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.interval = interval
        self.max_memory_mb = max_memory_mb
        self.max_cpu_percent = max_cpu_percent
        
        # Monitoring state
        self.monitoring = False
        self.monitor_thread = None
        self.start_time = None
        self.end_time = None
        
        # Data storage
        self.timestamps = []
        self.cpu_percent = []
        self.cpu_per_core = []
        self.memory_mb = []
        self.memory_percent = []
        self.disk_read_mb = []
        self.disk_write_mb = []
        self.net_sent_mb = []
        self.net_recv_mb = []
        
        # Get initial disk and network stats
        self.initial_disk_io = psutil.disk_io_counters()
        self.initial_net_io = psutil.net_io_counters()
        
        # Get process
        self.process = psutil.Process(os.getpid())
        
        logger.info(f"ResourceMonitor initialized: {output_dir}")
        logger.info(f"Sampling interval: {interval}s")
        if max_memory_mb:
            logger.info(f"Memory alert threshold: {max_memory_mb} MB")
        if max_cpu_percent:
            logger.info(f"CPU alert threshold: {max_cpu_percent}%")
    
    def start_monitoring(self):
        """Start monitoring resources in background thread."""
        if self.monitoring:
            logger.warning("Monitoring already active")
            return
        
        self.monitoring = True
        self.start_time = time.time()
        
        # Reset data
        self.timestamps = []
        self.cpu_percent = []
        self.cpu_per_core = []
        self.memory_mb = []
        self.memory_percent = []
        self.disk_read_mb = []
        self.disk_write_mb = []
        self.net_sent_mb = []
        self.net_recv_mb = []
        
        # Start monitoring thread
        self.monitor_thread = threading.Thread(target=self._monitor_loop, daemon=True)
        self.monitor_thread.start()
        
        logger.info("Started resource monitoring")
    
    def stop_monitoring(self) -> Dict:
        """
        Stop monitoring and return statistics.
        
        Returns
        -------
        dict
            Resource usage statistics
        """
        if not self.monitoring:
            logger.warning("Monitoring not active")
            return {}
        
        self.monitoring = False
        self.end_time = time.time()
        
        # Wait for thread to finish
        if self.monitor_thread:
            self.monitor_thread.join(timeout=5.0)
        
        # Calculate statistics
        stats = self._calculate_statistics()
        
        # Save statistics
        stats_file = self.output_dir / "resource_stats.json"
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2)
        
        logger.info("Stopped resource monitoring")
        logger.info(f"Statistics saved to: {stats_file}")
        
        return stats
    
    def _monitor_loop(self):
        """Main monitoring loop (runs in background thread)."""
        while self.monitoring:
            try:
                # Get timestamp
                timestamp = time.time() - self.start_time
                self.timestamps.append(timestamp)
                
                # CPU usage
                cpu_pct = psutil.cpu_percent(interval=None)
                self.cpu_percent.append(cpu_pct)
                
                cpu_per_core = psutil.cpu_percent(interval=None, percpu=True)
                self.cpu_per_core.append(cpu_per_core)
                
                # Memory usage
                mem = psutil.virtual_memory()
                mem_mb = mem.used / (1024 ** 2)
                self.memory_mb.append(mem_mb)
                self.memory_percent.append(mem.percent)
                
                # Disk I/O
                disk_io = psutil.disk_io_counters()
                if disk_io and self.initial_disk_io:
                    read_mb = (disk_io.read_bytes - self.initial_disk_io.read_bytes) / (1024 ** 2)
                    write_mb = (disk_io.write_bytes - self.initial_disk_io.write_bytes) / (1024 ** 2)
                    self.disk_read_mb.append(read_mb)
                    self.disk_write_mb.append(write_mb)
                
                # Network I/O
                net_io = psutil.net_io_counters()
                if net_io and self.initial_net_io:
                    sent_mb = (net_io.bytes_sent - self.initial_net_io.bytes_sent) / (1024 ** 2)
                    recv_mb = (net_io.bytes_recv - self.initial_net_io.bytes_recv) / (1024 ** 2)
                    self.net_sent_mb.append(sent_mb)
                    self.net_recv_mb.append(recv_mb)
                
                # Check thresholds and alert
                if self.max_memory_mb and mem_mb > self.max_memory_mb:
                    logger.warning(f"Memory usage ({mem_mb:.1f} MB) exceeded threshold ({self.max_memory_mb} MB)")
                
                if self.max_cpu_percent and cpu_pct > self.max_cpu_percent:
                    logger.warning(f"CPU usage ({cpu_pct:.1f}%) exceeded threshold ({self.max_cpu_percent}%)")
                
            except Exception as e:
                logger.error(f"Error in monitoring loop: {e}")
            
            # Sleep until next sample
            time.sleep(self.interval)
    
    def _calculate_statistics(self) -> Dict:
        """Calculate summary statistics from collected data."""
        import numpy as np
        
        runtime = self.end_time - self.start_time if self.end_time else 0
        
        stats = {
            'runtime': runtime,
            'n_samples': len(self.timestamps),
            'sampling_interval': self.interval,
            'cpu': {
                'mean': float(np.mean(self.cpu_percent)) if self.cpu_percent else 0,
                'median': float(np.median(self.cpu_percent)) if self.cpu_percent else 0,
                'std': float(np.std(self.cpu_percent)) if self.cpu_percent else 0,
                'min': float(np.min(self.cpu_percent)) if self.cpu_percent else 0,
                'max': float(np.max(self.cpu_percent)) if self.cpu_percent else 0,
                'n_cores': psutil.cpu_count()
            },
            'memory': {
                'mean_mb': float(np.mean(self.memory_mb)) if self.memory_mb else 0,
                'median_mb': float(np.median(self.memory_mb)) if self.memory_mb else 0,
                'std_mb': float(np.std(self.memory_mb)) if self.memory_mb else 0,
                'min_mb': float(np.min(self.memory_mb)) if self.memory_mb else 0,
                'peak_mb': float(np.max(self.memory_mb)) if self.memory_mb else 0,
                'total_mb': psutil.virtual_memory().total / (1024 ** 2)
            },
            'disk_io': {
                'total_read_mb': float(self.disk_read_mb[-1]) if self.disk_read_mb else 0,
                'total_write_mb': float(self.disk_write_mb[-1]) if self.disk_write_mb else 0,
                'read_rate_mbps': float(self.disk_read_mb[-1] / runtime) if self.disk_read_mb and runtime > 0 else 0,
                'write_rate_mbps': float(self.disk_write_mb[-1] / runtime) if self.disk_write_mb and runtime > 0 else 0
            },
            'network_io': {
                'total_sent_mb': float(self.net_sent_mb[-1]) if self.net_sent_mb else 0,
                'total_recv_mb': float(self.net_recv_mb[-1]) if self.net_recv_mb else 0,
                'send_rate_mbps': float(self.net_sent_mb[-1] / runtime) if self.net_sent_mb and runtime > 0 else 0,
                'recv_rate_mbps': float(self.net_recv_mb[-1] / runtime) if self.net_recv_mb and runtime > 0 else 0
            }
        }
        
        return stats
    
    def generate_plots(self):
        """Generate visualization plots for resource usage."""
        if not self.timestamps:
            logger.warning("No monitoring data available")
            return
        
        # Set style
        sns.set_style("whitegrid")
        
        # Create figure with subplots
        fig, axes = plt.subplots(2, 2, figsize=(15, 10))
        fig.suptitle('RAPTOR Resource Monitoring', fontsize=16, fontweight='bold')
        
        # Convert timestamps to minutes
        time_minutes = [t / 60 for t in self.timestamps]
        
        # Plot 1: CPU Usage
        axes[0, 0].plot(time_minutes, self.cpu_percent, color='#2E86AB', linewidth=2)
        axes[0, 0].axhline(y=np.mean(self.cpu_percent), color='r', linestyle='--', 
                          label=f'Mean: {np.mean(self.cpu_percent):.1f}%')
        axes[0, 0].fill_between(time_minutes, 0, self.cpu_percent, alpha=0.3, color='#2E86AB')
        axes[0, 0].set_xlabel('Time (minutes)', fontsize=12)
        axes[0, 0].set_ylabel('CPU Usage (%)', fontsize=12)
        axes[0, 0].set_title('CPU Usage Over Time', fontsize=14, fontweight='bold')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # Plot 2: Memory Usage
        axes[0, 1].plot(time_minutes, self.memory_mb, color='#A23B72', linewidth=2)
        axes[0, 1].axhline(y=np.mean(self.memory_mb), color='r', linestyle='--',
                          label=f'Mean: {np.mean(self.memory_mb):.0f} MB')
        axes[0, 1].fill_between(time_minutes, 0, self.memory_mb, alpha=0.3, color='#A23B72')
        axes[0, 1].set_xlabel('Time (minutes)', fontsize=12)
        axes[0, 1].set_ylabel('Memory Usage (MB)', fontsize=12)
        axes[0, 1].set_title('Memory Usage Over Time', fontsize=14, fontweight='bold')
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
        
        # Plot 3: Disk I/O
        if self.disk_read_mb and self.disk_write_mb:
            axes[1, 0].plot(time_minutes, self.disk_read_mb, color='#18A558', 
                           linewidth=2, label='Read')
            axes[1, 0].plot(time_minutes, self.disk_write_mb, color='#F77F00', 
                           linewidth=2, label='Write')
            axes[1, 0].set_xlabel('Time (minutes)', fontsize=12)
            axes[1, 0].set_ylabel('Data Transfer (MB)', fontsize=12)
            axes[1, 0].set_title('Disk I/O Over Time', fontsize=14, fontweight='bold')
            axes[1, 0].legend()
            axes[1, 0].grid(True, alpha=0.3)
        
        # Plot 4: Summary Statistics
        axes[1, 1].axis('off')
        stats = self._calculate_statistics()
        
        summary_text = f"""
        RESOURCE USAGE SUMMARY
        
        Runtime: {stats['runtime']/60:.1f} minutes
        Samples: {stats['n_samples']}
        
        CPU:
          Mean: {stats['cpu']['mean']:.1f}%
          Peak: {stats['cpu']['max']:.1f}%
          Cores: {stats['cpu']['n_cores']}
        
        Memory:
          Mean: {stats['memory']['mean_mb']:.0f} MB
          Peak: {stats['memory']['peak_mb']:.0f} MB
          Total: {stats['memory']['total_mb']:.0f} MB
        
        Disk I/O:
          Read: {stats['disk_io']['total_read_mb']:.1f} MB
          Write: {stats['disk_io']['total_write_mb']:.1f} MB
        
        Network:
          Sent: {stats['network_io']['total_sent_mb']:.1f} MB
          Received: {stats['network_io']['total_recv_mb']:.1f} MB
        """
        
        axes[1, 1].text(0.1, 0.5, summary_text, fontsize=11, verticalalignment='center',
                       family='monospace', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
        
        # Adjust layout and save
        plt.tight_layout()
        
        output_file = self.output_dir / "resource_usage.png"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        logger.info(f"Resource usage plot saved to: {output_file}")
        
        plt.close()
        
        # Generate individual plots
        self._generate_individual_plots(time_minutes)
    
    def _generate_individual_plots(self, time_minutes):
        """Generate individual plots for each metric."""
        # CPU plot
        plt.figure(figsize=(10, 6))
        plt.plot(time_minutes, self.cpu_percent, color='#2E86AB', linewidth=2)
        plt.xlabel('Time (minutes)', fontsize=12)
        plt.ylabel('CPU Usage (%)', fontsize=12)
        plt.title('CPU Usage', fontsize=14, fontweight='bold')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(self.output_dir / "cpu_usage.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # Memory plot
        plt.figure(figsize=(10, 6))
        plt.plot(time_minutes, self.memory_mb, color='#A23B72', linewidth=2)
        plt.xlabel('Time (minutes)', fontsize=12)
        plt.ylabel('Memory Usage (MB)', fontsize=12)
        plt.title('Memory Usage', fontsize=14, fontweight='bold')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(self.output_dir / "memory_usage.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Individual plots saved to: {self.output_dir}")
    
    def export_data(self) -> pd.DataFrame:
        """
        Export monitoring data as DataFrame.
        
        Returns
        -------
        pd.DataFrame
            Monitoring data
        """
        if not self.timestamps:
            logger.warning("No monitoring data available")
            return pd.DataFrame()
        
        data = {
            'timestamp': self.timestamps,
            'cpu_percent': self.cpu_percent,
            'memory_mb': self.memory_mb,
            'memory_percent': self.memory_percent
        }
        
        if self.disk_read_mb:
            data['disk_read_mb'] = self.disk_read_mb
            data['disk_write_mb'] = self.disk_write_mb
        
        if self.net_sent_mb:
            data['net_sent_mb'] = self.net_sent_mb
            data['net_recv_mb'] = self.net_recv_mb
        
        df = pd.DataFrame(data)
        
        # Save to CSV
        csv_file = self.output_dir / "monitoring_data.csv"
        df.to_csv(csv_file, index=False)
        logger.info(f"Monitoring data exported to: {csv_file}")
        
        return df


# =============================================================================
# Convenience Functions
# =============================================================================

def monitor_function(func, output_dir: str = "monitoring", **kwargs):
    """
    Decorator to monitor a function's resource usage.
    
    Parameters
    ----------
    func : callable
        Function to monitor
    output_dir : str
        Output directory for monitoring data
    **kwargs
        Additional arguments for ResourceMonitor
    
    Returns
    -------
    callable
        Decorated function
    
    Examples
    --------
    >>> @monitor_function
    ... def my_analysis():
    ...     # Your analysis code
    ...     pass
    """
    def wrapper(*args, **func_kwargs):
        monitor = ResourceMonitor(output_dir=output_dir, **kwargs)
        monitor.start_monitoring()
        
        try:
            result = func(*args, **func_kwargs)
            return result
        finally:
            stats = monitor.stop_monitoring()
            monitor.generate_plots()
            
            print(f"\n=== Resource Usage Summary ===")
            print(f"Runtime: {stats['runtime']:.1f}s")
            print(f"Peak CPU: {stats['cpu']['max']:.1f}%")
            print(f"Peak Memory: {stats['memory']['peak_mb']:.1f} MB")
            print(f"Results saved to: {output_dir}")
    
    return wrapper


if __name__ == '__main__':
    print("RAPTOR Resource Monitor")
    print("======================")
    print("\nMonitor CPU, memory, and I/O during pipeline execution.")
    print("\nUsage:")
    print("  from raptor.resource_monitoring import ResourceMonitor")
    print("  monitor = ResourceMonitor()")
    print("  monitor.start_monitoring()")
    print("  # ... run analysis ...")
    print("  stats = monitor.stop_monitoring()")
    print("  monitor.generate_plots()")
