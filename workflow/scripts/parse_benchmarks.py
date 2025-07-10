#!/usr/bin/env python3
"""
Extract key benchmark metrics from Snakemake benchmark files
"""

import os
import pandas as pd
import glob
from pathlib import Path

def parse_benchmark_file(filepath):
    """Parse a single benchmark file and extract key metrics"""
    try:
        # Read the benchmark file (tab-separated)
        df = pd.read_csv(filepath, sep='\t')
        
        # Get the first (and usually only) row
        row = df.iloc[0]
        
        # Extract key metrics
        metrics = {
            'file': os.path.basename(filepath),
            'sample': extract_sample_name(filepath),
            'tool': extract_tool_name(filepath),
            'wall_time_seconds': row['s'],
            'wall_time_formatted': row['h:m:s'],
            'peak_ram_gb': round(row['max_rss'] / 1024, 2),
            'peak_ram_mb': round(row['max_rss'], 2),
            'cpu_efficiency_percent': round(row['mean_load'], 2),
            'io_read_gb': round(row['io_in'] / 1024, 2),
            'io_write_gb': round(row['io_out'] / 1024, 2),
            'cpu_time_seconds': row['cpu_time']
        }
        
        return metrics
        
    except Exception as e:
        print(f"Error parsing {filepath}: {e}")
        return None

def extract_sample_name(filepath):
    """Extract sample name from file path"""
    # Assuming format: benchmark/SAMPLE/tool.txt
    path_parts = Path(filepath).parts
    if len(path_parts) >= 2:
        return path_parts[-2]  # Second to last part should be sample
    return "unknown"

def extract_tool_name(filepath):
    """Extract tool name from filename"""
    filename = os.path.basename(filepath)
    # Remove .txt extension and any sample prefixes
    tool_name = filename.replace('.txt', '')
    
    # Handle cases like "busco_048ss_spades.txt"
    if '_' in tool_name:
        parts = tool_name.split('_')
        # First part is usually the tool name
        return parts[0]
    
    return tool_name

def main():
    # Set the benchmark directory path
    benchmark_dir = "benchmark"  # Adjust path as needed
    
    if not os.path.exists(benchmark_dir):
        print(f"Benchmark directory '{benchmark_dir}' not found!")
        print("Please check the path or run from the correct directory.")
        return
    
    # Find all .txt files in benchmark directory
    pattern = os.path.join(benchmark_dir, "**", "*.txt")
    benchmark_files = glob.glob(pattern, recursive=True)
    
    if not benchmark_files:
        print(f"No benchmark files found in {benchmark_dir}")
        return
    
    print(f"Found {len(benchmark_files)} benchmark files")
    
    # Parse all benchmark files
    all_metrics = []
    for filepath in sorted(benchmark_files):
        metrics = parse_benchmark_file(filepath)
        if metrics:
            all_metrics.append(metrics)
    
    if not all_metrics:
        print("No valid benchmark data found!")
        return
    
    # Convert to DataFrame for easy manipulation
    df = pd.DataFrame(all_metrics)
    
    # Display summary
    print("\n" + "="*80)
    print("BENCHMARK SUMMARY")
    print("="*80)
    
    # Group by tool
    print("\nBy Tool:")
    print("-" * 40)
    tool_summary = df.groupby('tool').agg({
        'wall_time_seconds': ['count', 'mean', 'min', 'max'],
        'peak_ram_gb': ['mean', 'min', 'max'],
        'cpu_efficiency_percent': 'mean'
    }).round(2)
    print(tool_summary)
    
    # Show all individual results
    print("\n\nDetailed Results:")
    print("-" * 40)
    
    # Sort by tool then sample
    df_sorted = df.sort_values(['tool', 'sample'])
    
    for _, row in df_sorted.iterrows():
        print(f"{row['tool']:15} | {row['sample']:10} | "
              f"{row['wall_time_formatted']:>8} | "
              f"{row['peak_ram_gb']:>6.1f}GB | "
              f"{row['cpu_efficiency_percent']:>6.1f}% | "
              f"I/O: {row['io_read_gb']:>5.1f}R/{row['io_write_gb']:>5.1f}W GB")
    
    # Save detailed results to CSV
    output_file = "benchmark_summary.csv"
    df.to_csv(output_file, index=False)
    print(f"\n📊 Detailed results saved to: {output_file}")
    
    # Performance insights
    print("\n" + "="*80)
    print("PERFORMANCE INSIGHTS")
    print("="*80)
    
    # Fastest/slowest tools
    fastest = df.loc[df['wall_time_seconds'].idxmin()]
    slowest = df.loc[df['wall_time_seconds'].idxmax()]
    
    print(f"⚡ Fastest: {fastest['tool']} ({fastest['wall_time_formatted']})")
    print(f"🐌 Slowest: {slowest['tool']} ({slowest['wall_time_formatted']})")
    
    # Memory usage
    most_memory = df.loc[df['peak_ram_gb'].idxmax()]
    least_memory = df.loc[df['peak_ram_gb'].idxmin()]
    
    print(f"🧠 Most memory: {most_memory['tool']} ({most_memory['peak_ram_gb']} GB)")
    print(f"💾 Least memory: {least_memory['tool']} ({least_memory['peak_ram_gb']} GB)")
    
    # CPU efficiency
    most_efficient = df.loc[df['cpu_efficiency_percent'].idxmax()]
    least_efficient = df.loc[df['cpu_efficiency_percent'].idxmin()]
    
    print(f"🔥 Most CPU efficient: {most_efficient['tool']} ({most_efficient['cpu_efficiency_percent']}%)")
    print(f"😴 Least CPU efficient: {least_efficient['tool']} ({least_efficient['cpu_efficiency_percent']}%)")

if __name__ == "__main__":
    main()