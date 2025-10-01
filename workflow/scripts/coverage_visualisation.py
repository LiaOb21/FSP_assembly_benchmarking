#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse
import os

def parse_flagstat(flagstat_file):
    """
    Parse flagstat output to extract mapping statistics.
    
    Args:
        flagstat_file (str): Path to the flagstat output file
        
    Returns:
        dict: Dictionary containing mapping statistics
    """
    mapping_stats = {}
    
    try:
        with open(flagstat_file, 'r') as f:
            lines = f.readlines()
        
        # Parse each line looking for specific patterns
        for line in lines:
            line = line.strip()
            
            # Total reads: "13479 + 0 in total (QC-passed reads + QC-failed reads)"
            if 'in total' in line:
                mapping_stats['total_reads'] = int(line.split()[0])
                
            # Mapped reads: "2495 + 0 mapped (18.51% : N/A)" (but not "primary mapped")
            elif 'mapped (' in line and 'primary' not in line:
                parts = line.split()
                mapping_stats['mapped_reads'] = int(parts[0])
                # Extract percentage from format "mapped (18.51% : N/A)"
                pct_part = line.split('(')[1].split('%')[0]
                mapping_stats['mapping_rate'] = float(pct_part)
                
            # Properly paired reads: "2340 + 0 properly paired (17.38% : N/A)"
            elif 'properly paired' in line:
                mapping_stats['properly_paired'] = int(line.split()[0])
        
        # Calculate unmapped reads if we have the required data
        if 'total_reads' in mapping_stats and 'mapped_reads' in mapping_stats:
            mapping_stats['unmapped_reads'] = mapping_stats['total_reads'] - mapping_stats['mapped_reads']
            
    except Exception as e:
        print(f"Warning: Could not parse flagstat file: {e}")
        return {}
    
    return mapping_stats

def create_coverage_plots(coverage_file, output_file, sample_name, assembler_name, flagstat_file=None):
    """
    Create comprehensive coverage visualization plots with 6 subplots.
    
    Args:
        coverage_file (str): Path to samtools coverage output file
        output_file (str): Path for output PNG file
        sample_name (str): Name of the sample
        assembler_name (str): Name of the assembler
        flagstat_file (str, optional): Path to flagstat file for mapping statistics
        
    Returns:
        dict: Dictionary containing all calculated statistics
    """
    
    # ==================== DATA LOADING ====================
    
    # Read coverage data from samtools coverage output
    try:
        coverage_data = pd.read_csv(coverage_file, sep='\t')
        print(f"Loaded {len(coverage_data)} contigs from {coverage_file}")
    except Exception as e:
        print(f"Error reading file: {e}")
        return
    
    # Read mapping statistics if flagstat file is provided
    mapping_stats = {}
    if flagstat_file and os.path.exists(flagstat_file):
        mapping_stats = parse_flagstat(flagstat_file)
        if mapping_stats:  # Only print if we successfully parsed the file
            print(f"Loaded mapping statistics from {flagstat_file}")
        else:
            print(f"Warning: Could not parse mapping statistics from {flagstat_file}")
    
    # ==================== CALCULATE STATISTICS ====================
    
    # Basic coverage statistics
    avg_coverage = coverage_data['meandepth'].mean()
    median_coverage = coverage_data['meandepth'].median()
    peak_coverage = coverage_data['meandepth'].max()
    min_coverage = coverage_data['meandepth'].min()
    std_coverage = coverage_data['meandepth'].std()
    total_bases = coverage_data['endpos'].sum()
    
    # Weighted average coverage (larger contigs have more influence)
    weighted_avg_coverage = np.average(coverage_data['meandepth'], weights=coverage_data['endpos'])
    
    # Overall coverage percentage across all contigs
    covered_bases = (coverage_data['coverage'] * coverage_data['endpos'] / 100).sum()
    overall_coverage_pct = (covered_bases / total_bases) * 100 if total_bases > 0 else 0
    
    # Longest and shortest contigs
    longest_contig = coverage_data['endpos'].max()
    shortest_contig = coverage_data['endpos'].min()
    
    # Count contigs smaller than 250bp
    small_contigs_count = len(coverage_data[coverage_data['endpos'] < 250])
    small_contigs_pct = (small_contigs_count / len(coverage_data)) * 100

    # Count contigs longer than 1kb
    long_contigs_count = len(coverage_data[coverage_data['endpos'] >= 1000])
    long_contigs_pct = (long_contigs_count / len(coverage_data)) * 100

    # Count contigs between 250bp and 1kb
    mid_contigs_count = len(coverage_data[(coverage_data['endpos'] >= 250) & (coverage_data['endpos'] < 1000)])
    mid_contigs_pct = (mid_contigs_count / len(coverage_data)) * 100
    
    # ==================== PRINT SUMMARY STATISTICS ====================
    
    print(f"\n=== Coverage Summary for {sample_name} ({assembler_name}) ===")
    print(f"Total contigs: {len(coverage_data):,}")
    print(f"Contigs < 250bp: {small_contigs_count:,} ({small_contigs_pct:.1f}%)")
    print(f"Total bases: {total_bases:,}")
    print(f"Average coverage: {avg_coverage:.2f}x")
    print(f"Weighted average coverage: {weighted_avg_coverage:.2f}x")
    print(f"Median coverage: {median_coverage:.2f}x")
    print(f"Peak coverage: {peak_coverage:.2f}x")
    print(f"Min coverage: {min_coverage:.2f}x")
    print(f"Standard deviation: {std_coverage:.2f}x")
    print(f"Overall coverage: {overall_coverage_pct:.2f}%")
    
    # Print mapping statistics if available
    if mapping_stats:
        print(f"\n=== Mapping Statistics ===")
        print(f"Total reads: {mapping_stats.get('total_reads', 'N/A'):,}")
        print(f"Mapped reads: {mapping_stats.get('mapped_reads', 'N/A'):,}")
        print(f"Unmapped reads: {mapping_stats.get('unmapped_reads', 'N/A'):,}")
        print(f"Mapping rate: {mapping_stats.get('mapping_rate', 'N/A'):.1f}%")
        print(f"Properly paired: {mapping_stats.get('properly_paired', 'N/A'):,}")
    else:
        print(f"\n=== Mapping Statistics ===")
        print("No flagstat file provided - mapping statistics not available")
    
    # ==================== CREATE VISUALIZATION ====================
    
    # Create 2x3 subplot grid
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle(f'Coverage Analysis - {sample_name} ({assembler_name})', 
                fontsize=16, fontweight='bold')
    
    # Plot 1: Coverage distribution histogram (top-left)
    ax1 = axes[0, 0]
    n, bins, patches = ax1.hist(coverage_data['meandepth'], bins=50, alpha=0.7, 
                            color='skyblue', edgecolor='black', density=False)

    # Force Y-axis to show integer counts
    ax1.yaxis.set_major_locator(plt.MaxNLocator(integer=True))  # Force integer ticks
    ax1.set_ylim(0, max(n) + 1)  # Set Y-axis limit to actual max count + 1

    # Add vertical lines for key statistics
    ax1.axvline(avg_coverage, color='red', linestyle='--', linewidth=2, 
            label=f'Mean: {avg_coverage:.1f}x')
    ax1.axvline(median_coverage, color='green', linestyle='--', linewidth=2, 
            label=f'Median: {median_coverage:.1f}x')
    ax1.axvline(peak_coverage, color='orange', linestyle='--', linewidth=2, 
            label=f'Peak: {peak_coverage:.1f}x')
    ax1.set_xlabel('Coverage Depth (x)')
    ax1.set_ylabel('Number of Contigs')
    ax1.set_title('Coverage Distribution')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Coverage vs Contig Length scatter plot (top-middle)
    ax2 = axes[0, 1]
    # Color points by coverage percentage, size by default
    scatter = ax2.scatter(coverage_data['endpos'], coverage_data['meandepth'], 
                         alpha=0.6, c=coverage_data['coverage'], cmap='viridis', s=20)
    # Add horizontal reference lines
    ax2.axhline(avg_coverage, color='red', linestyle='--', linewidth=2, 
               label=f'Mean: {avg_coverage:.1f}x')
    ax2.axhline(weighted_avg_coverage, color='purple', linestyle='--', linewidth=2, 
               label=f'Weighted Mean: {weighted_avg_coverage:.1f}x')
    ax2.set_xlabel('Contig Length (bp)')
    ax2.set_ylabel('Coverage Depth (x)')
    ax2.set_title('Coverage vs Contig Length')
    ax2.set_xscale('log')  # Log scale for better visualization of length distribution
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    plt.colorbar(scatter, ax=ax2, label='Coverage %')
    
    # Plot 3: Cumulative coverage distribution (top-right)
    ax3 = axes[0, 2]
    # Sort coverage values and create cumulative count
    sorted_coverage = np.sort(coverage_data['meandepth'])
    cumulative_count = np.arange(1, len(sorted_coverage) + 1)
    ax3.plot(sorted_coverage, cumulative_count, linewidth=2, color='blue')
    # Add vertical reference lines
    ax3.axvline(avg_coverage, color='red', linestyle='--', linewidth=2, 
               label=f'Mean: {avg_coverage:.1f}x')
    ax3.axvline(median_coverage, color='green', linestyle='--', linewidth=2, 
               label=f'Median: {median_coverage:.1f}x')
    ax3.set_xlabel('Coverage Depth (x)')
    ax3.set_ylabel('Number of Contigs')
    ax3.set_title('Cumulative Coverage Distribution')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Contig length distribution (bottom-left)
    ax4 = axes[1, 0]
    ax4.hist(coverage_data['endpos'], bins=50, alpha=0.7, color='lightgreen', 
            edgecolor='black')
    # Add vertical lines for thresholds
    ax4.axvline(250, color='red', linestyle='--', linewidth=2, 
            label=f'250bp threshold\n({small_contigs_count} contigs below)')
    ax4.axvline(1000, color='blue', linestyle='--', linewidth=2, 
            label=f'1kb threshold\n({long_contigs_count} contigs above)')
    ax4.set_xlabel('Contig Length (bp)')
    ax4.set_ylabel('Number of Contigs')
    ax4.set_title('Contig Length Distribution')
    ax4.set_xscale('log')  # Log scale for better visualization
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    

    # Plot 5: Coverage Statistics text box (bottom-middle) - CENTERED VERSION
    ax5 = axes[1, 1]
    ax5.axis('off')  # Turn off axes for text display
    
    # Prepare comprehensive coverage statistics text
    coverage_stats_text = f"""Coverage Statistics:

Mean Coverage: {avg_coverage:.2f}x
Weighted Mean: {weighted_avg_coverage:.2f}x
Median Coverage: {median_coverage:.2f}x
Peak Coverage: {peak_coverage:.2f}x
Min Coverage: {min_coverage:.2f}x
Std Deviation: {std_coverage:.2f}x

Assembly Info:
Total Contigs: {len(coverage_data):,}
Contigs < 250bp: {small_contigs_count:,} ({small_contigs_pct:.1f}%)
Contigs 250bp-1kb: {mid_contigs_count:,} ({mid_contigs_pct:.1f}%)
Contigs > 1kb: {long_contigs_count:,} ({long_contigs_pct:.1f}%)
Total Length: {total_bases:,} bp
Longest Contig: {longest_contig:,} bp
Shortest Contig: {shortest_contig:,} bp
Overall Coverage: {overall_coverage_pct:.1f}%
"""
    
    # Display coverage statistics CENTERED in the quadrant
    ax5.text(0.5, 0.5, coverage_stats_text, transform=ax5.transAxes, 
            fontsize=12,  # Increased from 11
            verticalalignment='center',  # Changed from 'top' to 'center'
            horizontalalignment='center',  # Added horizontal centering
            fontfamily='monospace',
            bbox=dict(boxstyle='round,pad=0.8', facecolor='lightblue', alpha=0.9, 
                        edgecolor='darkblue', linewidth=1.5))  # Enhanced styling
    ax5.set_title('Coverage Statistics', fontweight='bold', fontsize=14)  # Increased title size
    
    # Plot 6: Quality & Mapping Metrics text box (bottom-right) - CENTERED VERSION
    ax6 = axes[1, 2]
    ax6.axis('off')  # Turn off axes for text display
    
    # Prepare quality and mapping metrics text (conditional on flagstat availability)
    if mapping_stats:
        quality_metrics_text = f"""Quality & Mapping Metrics:

Mapping Statistics:
Total Reads: {mapping_stats.get('total_reads', 'N/A'):,}
Mapped Reads: {mapping_stats.get('mapped_reads', 'N/A'):,}
Unmapped Reads: {mapping_stats.get('unmapped_reads', 'N/A'):,}
Mapping Rate: {mapping_stats.get('mapping_rate', 'N/A'):.1f}%
Properly Paired: {mapping_stats.get('properly_paired', 'N/A'):,}

Quality Scores:
Mean Base Quality: {coverage_data['meanbaseq'].mean():.1f}
Median Base Quality: {coverage_data['meanbaseq'].median():.1f}
Mean Mapping Quality: {coverage_data['meanmapq'].mean():.1f}
Median Mapping Quality: {coverage_data['meanmapq'].median():.1f}

Coverage Info:
Total Mapped Reads: {coverage_data['numreads'].sum():,}
Avg Reads per Contig: {coverage_data['numreads'].mean():.1f}
"""
    else:
        # Fallback text when no flagstat file is provided
        quality_metrics_text = f"""Quality Metrics:
(No flagstat file provided)

Mean Base Quality: {coverage_data['meanbaseq'].mean():.1f}
Median Base Quality: {coverage_data['meanbaseq'].median():.1f}
Mean Mapping Quality: {coverage_data['meanmapq'].mean():.1f}
Median Mapping Quality: {coverage_data['meanmapq'].median():.1f}

Per-Contig Ranges:
Base Quality: {coverage_data['meanbaseq'].min():.1f} - {coverage_data['meanbaseq'].max():.1f}
Mapping Quality: {coverage_data['meanmapq'].min():.1f} - {coverage_data['meanmapq'].max():.1f}

Read Mapping:
Total Mapped Reads: {coverage_data['numreads'].sum():,}
Avg Reads per Contig: {coverage_data['numreads'].mean():.1f}
Max Reads per Contig: {coverage_data['numreads'].max():,}
"""
    
    # Display quality metrics CENTERED in the quadrant
    ax6.text(0.5, 0.5, quality_metrics_text, transform=ax6.transAxes, 
            fontsize=12,  # Increased from 10
            verticalalignment='center',  # Changed from 'top' to 'center'
            horizontalalignment='center',  # Added horizontal centering
            fontfamily='monospace',
            bbox=dict(boxstyle='round,pad=0.8', facecolor='lightcoral', alpha=0.9,
                        edgecolor='darkred', linewidth=1.5))  # Enhanced styling
    ax6.set_title('Quality & Mapping Metrics', fontweight='bold', fontsize=14)  # Increased title size
    
    # ==================== FINALIZE AND SAVE ====================
    
    # Adjust layout to prevent overlapping
    plt.tight_layout()
    
    # Save the plot with high resolution
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()  # Close the figure to free memory
    
    print(f"Plot saved to: {output_file}")
    
    # ==================== RETURN STATISTICS DICTIONARY ====================
    
    # Prepare comprehensive statistics dictionary for summary file
    stats_dict = {
        'sample': sample_name,
        'assembler': assembler_name,
        'total_contigs': len(coverage_data),
        'contigs_under_250bp': small_contigs_count,
        'contigs_under_250bp_percent': small_contigs_pct,
        'contigs_250bp_to_1kb': mid_contigs_count,
        'contigs_250bp_to_1kb_percent': mid_contigs_pct,
        'contigs_over_1kb': long_contigs_count,
        'contigs_over_1kb_percent': long_contigs_pct,
        'total_bases': total_bases,
        'mean_coverage': avg_coverage,
        'weighted_mean_coverage': weighted_avg_coverage,
        'median_coverage': median_coverage,
        'peak_coverage': peak_coverage,
        'min_coverage': min_coverage,
        'std_coverage': std_coverage,
        'overall_coverage_pct': overall_coverage_pct,
        'longest_contig': longest_contig,
        'shortest_contig': shortest_contig,
        'mean_base_quality': coverage_data['meanbaseq'].mean(),
        'median_base_quality': coverage_data['meanbaseq'].median(),
        'mean_mapping_quality': coverage_data['meanmapq'].mean(),
        'median_mapping_quality': coverage_data['meanmapq'].median(),
    }
    
    # Add mapping statistics if available
    if mapping_stats:
        stats_dict.update({
            'total_reads': mapping_stats.get('total_reads', 'N/A'),
            'mapped_reads': mapping_stats.get('mapped_reads', 'N/A'),
            'unmapped_reads': mapping_stats.get('unmapped_reads', 'N/A'),
            'mapping_rate': mapping_stats.get('mapping_rate', 'N/A'),
            'properly_paired': mapping_stats.get('properly_paired', 'N/A')
        })
    
    return stats_dict

def main():
    parser = argparse.ArgumentParser(description='Create coverage visualization plots')
    parser.add_argument('coverage_file', help='Input coverage file (samtools coverage output)')
    parser.add_argument('output_file', help='Output plot file (PNG/PDF)')
    parser.add_argument('--sample', default='Sample', help='Sample name')
    parser.add_argument('--assembler', default='Assembler', help='Assembler name')
    parser.add_argument('--flagstat', help='Optional flagstat file for mapping statistics')
    parser.add_argument('--summary', help='Output summary file path')  # Add this line
    
    args = parser.parse_args()
    
    # Create the visualization
    stats = create_coverage_plots(args.coverage_file, args.output_file, 
                                 args.sample, args.assembler, args.flagstat)
    
    # Use provided summary file path or generate one
    if args.summary:
        summary_file = args.summary
    else:
        summary_file = args.output_file.replace('.png', '_summary.txt').replace('.pdf', '_summary.txt')
    
    # Save summary file
    with open(summary_file, 'w') as f:
        f.write(f"Coverage Analysis Summary - {args.sample} ({args.assembler})\n")
        f.write("=" * 60 + "\n\n")
        for key, value in stats.items():
            if isinstance(value, float):
                f.write(f"{key}: {value:.2f}\n")
            else:
                f.write(f"{key}: {value}\n")
    
    print(f"Summary statistics saved to: {summary_file}")

# Run the main function if script is executed directly
if __name__ == "__main__":
    main()