#!/usr/bin/env python3
"""
RAPTOR v2.1.1 Example Script: Resource Monitoring

Demonstrates real-time resource tracking during pipeline execution:
- CPU usage monitoring (overall and per-core)
- Memory consumption tracking
- Disk I/O monitoring
- Network activity tracking
- Bottleneck detection
- Performance visualization

Author: Ayeh Bolouki
License: MIT
"""

import argparse
import json
import sys
import time
import threading
from datetime import datetime
from pathlib import Path

# Check for dependencies
try:
    import numpy as np
except ImportError:
    print("ERROR: numpy is required")
    print("Install with: pip install numpy")
    sys.exit(1)

# Check for psutil
PSUTIL_AVAILABLE = True
try:
    import psutil
except ImportError:
    PSUTIL_AVAILABLE = False

# RAPTOR imports with fallback
RAPTOR_AVAILABLE = True
try:
    from raptor.resource_monitoring import ResourceMonitor, monitor_function
except ImportError:
    RAPTOR_AVAILABLE = False
    print("NOTE: RAPTOR modules not available. Running in demo mode.")


def print_banner():
    """Print RAPTOR banner."""
    print("""
    ╔══════════════════════════════════════════════════════════════╗
    ║       🦖 RAPTOR v2.1.1 - Resource Monitoring                ║
    ║                                                              ║
    ║   Real-Time CPU/Memory/Disk Tracking During Analysis        ║
    ╚══════════════════════════════════════════════════════════════╝
    """)


class SimpleResourceMonitor:
    """Simple resource monitor for demo mode or when RAPTOR not available."""
    
    def __init__(self, interval=1.0):
        self.interval = interval
        self.monitoring = False
        self.monitor_thread = None
        self.start_time = None
        self.end_time = None
        
        # Data storage
        self.timestamps = []
        self.cpu_percent = []
        self.memory_mb = []
        self.memory_percent = []
    
    def start_monitoring(self):
        """Start monitoring in background thread."""
        self.monitoring = True
        self.start_time = time.time()
        
        # Reset data
        self.timestamps = []
        self.cpu_percent = []
        self.memory_mb = []
        self.memory_percent = []
        
        # Start thread
        self.monitor_thread = threading.Thread(target=self._monitor_loop, daemon=True)
        self.monitor_thread.start()
        
        print(f"📊 Started resource monitoring (interval: {self.interval}s)")
    
    def stop_monitoring(self):
        """Stop monitoring and return statistics."""
        self.monitoring = False
        self.end_time = time.time()
        
        if self.monitor_thread:
            self.monitor_thread.join(timeout=5.0)
        
        return self._calculate_statistics()
    
    def _monitor_loop(self):
        """Main monitoring loop."""
        while self.monitoring:
            try:
                timestamp = time.time() - self.start_time
                self.timestamps.append(timestamp)
                
                if PSUTIL_AVAILABLE:
                    # Real monitoring
                    cpu = psutil.cpu_percent(interval=None)
                    mem = psutil.virtual_memory()
                    
                    self.cpu_percent.append(cpu)
                    self.memory_mb.append(mem.used / (1024 ** 2))
                    self.memory_percent.append(mem.percent)
                else:
                    # Simulated monitoring
                    self.cpu_percent.append(30 + np.random.normal(0, 15))
                    self.memory_mb.append(4000 + np.random.normal(0, 500))
                    self.memory_percent.append(40 + np.random.normal(0, 10))
                
            except Exception as e:
                print(f"Warning: Monitoring error: {e}")
            
            time.sleep(self.interval)
    
    def _calculate_statistics(self):
        """Calculate summary statistics."""
        runtime = self.end_time - self.start_time if self.end_time else 0
        
        if not self.cpu_percent:
            return {'runtime': runtime, 'n_samples': 0}
        
        stats = {
            'runtime': runtime,
            'n_samples': len(self.timestamps),
            'sampling_interval': self.interval,
            'cpu': {
                'mean': float(np.mean(self.cpu_percent)),
                'median': float(np.median(self.cpu_percent)),
                'std': float(np.std(self.cpu_percent)),
                'min': float(np.min(self.cpu_percent)),
                'max': float(np.max(self.cpu_percent)),
                'n_cores': psutil.cpu_count() if PSUTIL_AVAILABLE else 4
            },
            'memory': {
                'mean_mb': float(np.mean(self.memory_mb)),
                'median_mb': float(np.median(self.memory_mb)),
                'std_mb': float(np.std(self.memory_mb)),
                'min_mb': float(np.min(self.memory_mb)),
                'peak_mb': float(np.max(self.memory_mb)),
                'mean_percent': float(np.mean(self.memory_percent)),
                'total_mb': psutil.virtual_memory().total / (1024 ** 2) if PSUTIL_AVAILABLE else 16000
            }
        }
        
        return stats


def generate_demo_metrics(duration=30):
    """Generate demonstration metrics with simulated workload."""
    np.random.seed(42)
    
    n_samples = int(duration)
    timestamps = list(range(n_samples))
    
    # Simulate a pipeline with phases
    cpu_percent = []
    memory_mb = []
    
    for t in range(n_samples):
        # Phase 1: Data loading (low CPU, increasing memory)
        if t < n_samples * 0.1:
            cpu = 20 + np.random.normal(0, 5)
            mem = 2000 + t * 100 + np.random.normal(0, 100)
        
        # Phase 2: Alignment (high CPU, stable memory)
        elif t < n_samples * 0.4:
            cpu = 85 + np.random.normal(0, 10)
            mem = 4000 + np.random.normal(0, 200)
        
        # Phase 3: Quantification (medium CPU, moderate memory)
        elif t < n_samples * 0.6:
            cpu = 60 + np.random.normal(0, 15)
            mem = 6000 + np.random.normal(0, 300)
        
        # Phase 4: Statistics (low CPU, high memory)
        elif t < n_samples * 0.8:
            cpu = 30 + np.random.normal(0, 10)
            mem = 8000 + np.random.normal(0, 400)
        
        # Phase 5: Writing results (low CPU, decreasing memory)
        else:
            cpu = 15 + np.random.normal(0, 5)
            mem = 8000 - (t - n_samples * 0.8) * 200 + np.random.normal(0, 200)
        
        cpu_percent.append(max(0, min(100, cpu)))
        memory_mb.append(max(1000, mem))
    
    return {
        'runtime': duration,
        'n_samples': n_samples,
        'sampling_interval': 1.0,
        'timestamps': timestamps,
        'cpu_percent': cpu_percent,
        'memory_mb': memory_mb,
        'cpu': {
            'mean': float(np.mean(cpu_percent)),
            'median': float(np.median(cpu_percent)),
            'std': float(np.std(cpu_percent)),
            'min': float(np.min(cpu_percent)),
            'max': float(np.max(cpu_percent)),
            'n_cores': 8
        },
        'memory': {
            'mean_mb': float(np.mean(memory_mb)),
            'median_mb': float(np.median(memory_mb)),
            'std_mb': float(np.std(memory_mb)),
            'min_mb': float(np.min(memory_mb)),
            'peak_mb': float(np.max(memory_mb)),
            'mean_percent': float(np.mean(memory_mb) / 16000 * 100),
            'total_mb': 16000.0
        },
        'disk_io': {
            'total_read_mb': 2500.0,
            'total_write_mb': 850.0
        },
        'network_io': {
            'total_sent_mb': 15.0,
            'total_recv_mb': 45.0
        }
    }


def display_resource_stats(stats):
    """Display resource statistics with formatting."""
    
    print("\n" + "="*70)
    print("  🦖 RESOURCE MONITORING SUMMARY")
    print("="*70)
    
    # Runtime
    runtime = stats['runtime']
    print(f"\n  ⏱️  Runtime: {runtime:.1f} seconds ({runtime/60:.1f} minutes)")
    print(f"      Samples collected: {stats['n_samples']}")
    
    # CPU Statistics
    cpu = stats['cpu']
    print(f"\n  💻 CPU Usage:")
    
    # CPU bar
    cpu_mean = cpu['mean']
    bar_len = int(cpu_mean / 100 * 40)
    bar = '█' * bar_len + '░' * (40 - bar_len)
    print(f"      Mean:   [{bar}] {cpu_mean:5.1f}%")
    
    bar_len = int(cpu['max'] / 100 * 40)
    bar = '█' * bar_len + '░' * (40 - bar_len)
    print(f"      Peak:   [{bar}] {cpu['max']:5.1f}%")
    
    print(f"      Std:    {cpu['std']:.1f}%")
    print(f"      Cores:  {cpu['n_cores']}")
    
    # Memory Statistics
    mem = stats['memory']
    print(f"\n  🧠 Memory Usage:")
    
    # Memory bar (as percentage of total)
    mem_pct = mem['mean_mb'] / mem['total_mb'] * 100
    bar_len = int(mem_pct / 100 * 40)
    bar = '█' * bar_len + '░' * (40 - bar_len)
    print(f"      Mean:   [{bar}] {mem['mean_mb']/1024:.1f} GB ({mem_pct:.1f}%)")
    
    peak_pct = mem['peak_mb'] / mem['total_mb'] * 100
    bar_len = int(peak_pct / 100 * 40)
    bar = '█' * bar_len + '░' * (40 - bar_len)
    print(f"      Peak:   [{bar}] {mem['peak_mb']/1024:.1f} GB ({peak_pct:.1f}%)")
    
    print(f"      Total:  {mem['total_mb']/1024:.1f} GB")
    
    # Disk I/O (if available)
    if 'disk_io' in stats:
        disk = stats['disk_io']
        print(f"\n  💾 Disk I/O:")
        print(f"      Read:   {disk['total_read_mb']:.1f} MB")
        print(f"      Write:  {disk['total_write_mb']:.1f} MB")
    
    # Network I/O (if available)
    if 'network_io' in stats:
        net = stats['network_io']
        print(f"\n  🌐 Network I/O:")
        print(f"      Sent:   {net['total_sent_mb']:.1f} MB")
        print(f"      Recv:   {net['total_recv_mb']:.1f} MB")
    
    # Recommendations
    print(f"\n  📋 Recommendations:")
    
    if cpu['max'] > 90:
        print(f"      ⚠️  High CPU usage detected. Consider reducing parallelism.")
    elif cpu['mean'] < 30:
        print(f"      💡 Low CPU utilization. Consider increasing parallelism.")
    else:
        print(f"      ✓ CPU usage within normal range.")
    
    if peak_pct > 80:
        print(f"      ⚠️  High memory usage ({peak_pct:.1f}%). Monitor for OOM issues.")
    elif peak_pct < 30:
        print(f"      💡 Low memory usage. Could handle larger datasets.")
    else:
        print(f"      ✓ Memory usage within normal range.")


def run_resource_monitoring(command=None, analyze_file=None, interval=1.0, 
                            duration=30, demo=False):
    """Run resource monitoring."""
    
    if demo or not RAPTOR_AVAILABLE:
        print("\n🎮 Running in DEMO mode with simulated workload...")
        print(f"   Simulating {duration} seconds of pipeline execution...")
        
        # Generate demo metrics
        stats = generate_demo_metrics(duration=duration)
        
    elif analyze_file:
        # Analyze existing resource log
        print(f"\n📂 Loading resource log: {analyze_file}")
        with open(analyze_file, 'r') as f:
            stats = json.load(f)
        
    elif command:
        # Monitor a command
        print(f"\n🔍 Monitoring command: {command}")
        
        if RAPTOR_AVAILABLE:
            monitor = ResourceMonitor(output_dir='resource_monitoring', interval=interval)
            monitor.start_monitoring()
            
            # Execute command
            import subprocess
            result = subprocess.run(command, shell=True)
            
            stats = monitor.stop_monitoring()
            monitor.generate_plots()
        else:
            # Use simple monitor
            monitor = SimpleResourceMonitor(interval=interval)
            monitor.start_monitoring()
            
            import subprocess
            result = subprocess.run(command, shell=True)
            
            stats = monitor.stop_monitoring()
    
    else:
        # Interactive monitoring for specified duration
        print(f"\n📊 Starting {duration}-second monitoring session...")
        
        if RAPTOR_AVAILABLE:
            monitor = ResourceMonitor(output_dir='resource_monitoring', interval=interval)
        else:
            monitor = SimpleResourceMonitor(interval=interval)
        
        monitor.start_monitoring()
        
        # Show progress
        for i in range(duration):
            progress = (i + 1) / duration
            bar_len = int(progress * 40)
            bar = '█' * bar_len + '░' * (40 - bar_len)
            print(f"\r   Monitoring: [{bar}] {i+1}/{duration}s", end='', flush=True)
            time.sleep(1)
        
        print()  # New line after progress bar
        
        stats = monitor.stop_monitoring()
        
        if RAPTOR_AVAILABLE:
            monitor.generate_plots()
    
    # Display results
    display_resource_stats(stats)
    
    # Prepare output
    output = {
        'timestamp': datetime.now().isoformat(),
        'raptor_version': '2.1.1',
        'monitoring_config': {
            'interval': interval,
            'duration': duration if not command else stats.get('runtime', 0)
        },
        'statistics': stats
    }
    
    return output


def main():
    parser = argparse.ArgumentParser(
        description='🦖 RAPTOR Resource Monitoring',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Monitor a command
  python 09_resource_monitoring.py --monitor-command "python pipeline.py"
  
  # Monitor for a specific duration
  python 09_resource_monitoring.py --duration 60 --interval 0.5
  
  # Analyze existing resource log
  python 09_resource_monitoring.py --analyze resource_log.json
  
  # Demo mode (no actual monitoring)
  python 09_resource_monitoring.py --demo --duration 30

Features:
  - Real-time CPU/memory tracking (<1% overhead)
  - Per-process breakdown
  - Bottleneck detection
  - Resource prediction
  - Visualization generation
        """
    )
    
    parser.add_argument('--monitor-command', help='Command to monitor')
    parser.add_argument('--analyze', help='Analyze existing resource log')
    parser.add_argument('--interval', type=float, default=1.0, 
                        help='Sampling interval in seconds')
    parser.add_argument('--duration', type=int, default=30,
                        help='Monitoring duration in seconds')
    parser.add_argument('--output', '-o', help='Output JSON file')
    parser.add_argument('--demo', action='store_true', help='Run in demo mode')
    
    args = parser.parse_args()
    
    print_banner()
    
    # Check psutil availability
    if not PSUTIL_AVAILABLE and not args.demo and not args.analyze:
        print("⚠️  psutil not available - using simulated data")
        print("   Install with: pip install psutil")
    
    # Run monitoring
    results = run_resource_monitoring(
        command=args.monitor_command,
        analyze_file=args.analyze,
        interval=args.interval,
        duration=args.duration,
        demo=args.demo
    )
    
    # Save output
    output_file = args.output or 'resource_metrics.json'
    
    # Clean up non-serializable items
    clean_results = json.loads(json.dumps(results, default=lambda x: str(x) if not isinstance(x, (int, float, str, bool, list, dict, type(None))) else x))
    
    with open(output_file, 'w') as f:
        json.dump(clean_results, f, indent=2)
    
    print(f"\n💾 Results saved to: {output_file}")
    print("\n" + "="*70)
    print("  Making free science for everybody around the world 🌍")
    print("="*70 + "\n")


if __name__ == '__main__':
    main()
