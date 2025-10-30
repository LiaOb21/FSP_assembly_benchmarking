#!/usr/bin/env python3

import pandas as pd
import re
import os
from pathlib import Path
import glob
import argparse
import sys

def extract_sample_info(busco_filename, best_assembly_source_path):
    """
    Extract sample ID, assembler, reads type, kmer strategy and lineage info from multiple sources.
    
    Args:
        busco_filename (str): BUSCO summary filename
        best_assembly_source_path (str): Path to best_assembly_source.txt file
        
    Returns:
        tuple: (sample_id, lineage, reads_type, kmer_strategy, assembler)
    """
    # Extract sample_id and lineage from BUSCO filename
    # short_summary.specific.fungi_odb12.BUSCO_048ds_best_assembly_pilon.fa.txt
    pattern1 = r'short_summary\.specific\.([^.]+)\.BUSCO_(.+)_best_assembly_pilon\.fa\.txt'
    match1 = re.search(pattern1, busco_filename)
    
    if match1:
        lineage = match1.group(1)
        sample_id = match1.group(2)
        
        # Extract reads_type, kmer_strategy, and assembler from best_assembly_source.txt
        reads_type = None
        kmer_strategy = None
        assembler = None
        
        try:
            with open(best_assembly_source_path, 'r') as f:
                content = f.read().strip()
                # Format: "06_EG_R1R2_reads_length_2samples_SP_Lia:masurca"
                # or:     "06_EG_merged_kmergenie_2samples_SP_Lia:masurca"
                # or:     "06_EG_R1R2_manual_2samples_SP_Lia:masurca"
                
                # Split by colon to separate run info from assembler
                if ':' in content:
                    run_info, assembler = content.split(':', 1)
                    
                    # Extract reads type from run info
                    # Look for either R1R2 or merged
                    if 'R1R2' in run_info:
                        reads_type = 'R1R2'
                    elif 'merged' in run_info:
                        reads_type = 'merged'
                    
                    # Extract kmer strategy from run info
                    # Look for reads_length, kmergenie, or manual
                    if 'reads_length' in run_info:
                        kmer_strategy = 'reads_length'
                    elif 'kmergenie' in run_info:
                        kmer_strategy = 'kmergenie'
                    elif 'manual' in run_info:
                        kmer_strategy = 'manual'
                
        except FileNotFoundError:
            print(f"Warning: best_assembly_source.txt not found at {best_assembly_source_path}")
        except Exception as e:
            print(f"Error reading best_assembly_source.txt: {e}")
        
        return sample_id, lineage, reads_type, kmer_strategy, assembler

    return None, None, None, None, None


def parse_busco_file(filepath):
    """
    Parse a BUSCO summary file and extract relevant statistics.
    Excludes assembly statistics (num_scaffolds, num_contigs, total_length).
    
    Args:
        filepath (str): Path to BUSCO summary file
        
    Returns:
        dict: Dictionary containing parsed BUSCO statistics
    """
    results = {}
    
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Extract the results line (contains percentages)
    results_pattern = r'C:(\d+\.?\d*)%\[S:(\d+\.?\d*)%,D:(\d+\.?\d*)%\],F:(\d+\.?\d*)%,M:(\d+\.?\d*)%,n:(\d+)'
    match = re.search(results_pattern, content)
    
    if match:
        results['complete_percent'] = float(match.group(1))
        results['single_copy_percent'] = float(match.group(2))
        results['duplicated_percent'] = float(match.group(3))
        results['fragmented_percent'] = float(match.group(4))
        results['missing_percent'] = float(match.group(5))
        results['total_buscos'] = int(match.group(6))
    
    # Extract absolute numbers
    complete_match = re.search(r'(\d+)\s+Complete BUSCOs \(C\)', content)
    single_copy_match = re.search(r'(\d+)\s+Complete and single-copy BUSCOs \(S\)', content)
    duplicated_match = re.search(r'(\d+)\s+Complete and duplicated BUSCOs \(D\)', content)
    fragmented_match = re.search(r'(\d+)\s+Fragmented BUSCOs \(F\)', content)
    missing_match = re.search(r'(\d+)\s+Missing BUSCOs \(M\)', content)
    
    if complete_match:
        results['complete_count'] = int(complete_match.group(1))
    if single_copy_match:
        results['single_copy_count'] = int(single_copy_match.group(1))
    if duplicated_match:
        results['duplicated_count'] = int(duplicated_match.group(1))
    if fragmented_match:
        results['fragmented_count'] = int(fragmented_match.group(1))
    if missing_match:
        results['missing_count'] = int(missing_match.group(1))
    
    # REMOVED: Assembly statistics extraction
    # We no longer extract num_scaffolds, num_contigs, total_length
    
    return results


def parse_quast_file(filepath):
    """
    Parse a QUAST transposed_report.tsv file and extract relevant statistics.
    """
    try:
        df = pd.read_csv(filepath, sep='\t')
        if len(df) == 0:
            return {}
        
        row = df.iloc[0]
        
        results = {
            'total_contigs>250bp': row.get('# contigs', None),
            'contigs>1000bp': row.get('# contigs (>= 1000 bp)', None),
            'contigs>5000bp': row.get('# contigs (>= 5000 bp)', None),
            'contigs>10000bp': row.get('# contigs (>= 10000 bp)', None),
            'contigs>25000bp': row.get('# contigs (>= 25000 bp)', None),
            'total_length>250bp': row.get('Total length', None),
            'largest_contig': row.get('Largest contig', None),
            'gc_percent': row.get('GC (%)', None),
            'aUN': row.get('auN', None),
            'N50': row.get('N50', None),
            'N90': row.get('N90', None),
            'L50': row.get('L50', None),
            'L90': row.get('L90', None),
            'Ns_per_100kbp': row.get('Ns per 100 kbp', None)
        }

        # Convert to appropriate types - preserve decimals
        for key, value in results.items():
            if pd.isna(value):
                results[key] = None
            else:
                try:
                    # Try float first to preserve decimals
                    float_val = float(value)
                    # Only convert to int if it's actually a whole number
                    if float_val.is_integer():
                        results[key] = int(float_val)
                    else:
                        results[key] = float_val
                except (ValueError, TypeError):
                    results[key] = value
        
        return results
    except Exception as e:
        print(f"Error parsing QUAST file {filepath}: {e}")
        return {}
    

def parse_merqury_qv_file(filepath):
    """
    Parse a merquryfk.qv file and extract Error % and QV.
    
    Args:
        filepath (str): Path to merquryfk.qv file
        
    Returns:
        dict: Dictionary containing Error % and QV
    """
    try:
        df = pd.read_csv(filepath, sep='\t')
        if len(df) == 0:
            return {}
        
        row = df.iloc[0]  # Take first row
        
        results = {
            'error_percent': row.get('Error %', None),
            'qv': row.get('QV', None)
        }
        
        # Convert to appropriate types - preserve decimals
        for key, value in results.items():
            if pd.isna(value) or value == 'inf':
                results[key] = None if value != 'inf' else float('inf')
            else:
                try:
                    float_val = float(value)
                    results[key] = float_val
                except (ValueError, TypeError):
                    results[key] = value
        
        return results
        
    except Exception as e:
        print(f"Error parsing Merqury QV file {filepath}: {e}")
        return {}

def parse_merqury_completeness_file(filepath):
    """
    Parse a merquryfk.completeness.stats file and extract % Covered.
    
    Args:
        filepath (str): Path to merquryfk.completeness.stats file
        
    Returns:
        dict: Dictionary containing % Covered
    """
    try:
        df = pd.read_csv(filepath, sep='\t')
        if len(df) == 0:
            return {}
        
        row = df.iloc[0]  # Take first row
        
        results = {
            'percent_covered': row.get('% Covered', None)
        }
        
        # Convert to appropriate types - preserve decimals
        for key, value in results.items():
            if pd.isna(value):
                results[key] = None
            else:
                try:
                    float_val = float(value)
                    results[key] = float_val
                except (ValueError, TypeError):
                    results[key] = value
        
        return results
        
    except Exception as e:
        print(f"Error parsing Merqury completeness file {filepath}: {e}")
        return {}
    

def parse_coverage_file(filepath):
    """
    Parse a coverage summary file and extract relevant statistics.
    
    Args:
        filepath (str): Path to coverage summary file
        
    Returns:
        dict: Dictionary containing coverage statistics
    """
    try:
        results = {}
        
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if ':' in line:
                    key, value = line.split(':', 1)
                    key = key.strip()
                    value = value.strip()
                    
                    # Extract only the metrics we want
                    if key in ['mean_coverage', 'weighted_mean_coverage', 'median_coverage', 
                              'peak_coverage', 'min_coverage', 'std_coverage',
                              'mean_mapping_quality', 'median_mapping_quality',
                              'total_reads', 'mapped_reads', 'unmapped_reads',
                              'mapping_rate', 'properly_paired']:
                        
                        # Convert to appropriate types - preserve decimals
                        try:
                            float_val = float(value)
                            # Only convert to int if it's actually a whole number
                            if float_val.is_integer():
                                results[key] = int(float_val)
                            else:
                                results[key] = float_val
                        except (ValueError, TypeError):
                            results[key] = value
        
        return results
        
    except Exception as e:
        print(f"Error parsing coverage file {filepath}: {e}")
        return {}
    

def collect_all_metrics(base_dir, verbose=True):
    """
    Collection of BUSCO, QUAST, Merqury, and Coverage data.
    Collects both general and specific BUSCO results plus all other metrics.
    """
    data = []
    
    # Get all sample directories
    sample_dirs = [d for d in os.listdir(base_dir) 
                   if os.path.isdir(os.path.join(base_dir, d))]
    
    for sample_dir_name in sample_dirs:
        sample_path = os.path.join(base_dir, sample_dir_name)
        
        # Get QUAST data (should be same for both general and specific)
        quast_file = os.path.join(sample_path, 'best_assembly_info_and_QC', 
                                 'quast_pilon', 'transposed_report.tsv')
        quast_data = parse_quast_file(quast_file) if os.path.exists(quast_file) else {}
        
        # Get Merqury data (should be same for both general and specific)
        merqury_qv_file = os.path.join(sample_path, 'best_assembly_info_and_QC', 
                                      'merquryfk_pilon', 'merquryfk.qv')
        merqury_completeness_file = os.path.join(sample_path, 'best_assembly_info_and_QC', 
                                               'merquryfk_pilon', 'merquryfk.completeness.stats')
        
        merqury_qv_data = parse_merqury_qv_file(merqury_qv_file) if os.path.exists(merqury_qv_file) else {}
        merqury_completeness_data = parse_merqury_completeness_file(merqury_completeness_file) if os.path.exists(merqury_completeness_file) else {}
        
        # Get Coverage data (should be same for both general and specific)
        coverage_file = os.path.join(sample_path, 'best_assembly_info_and_QC', 
                                    'coverage_viz_pilon', 
                                    f'{sample_dir_name}_best_assembly_pilon_coverage_summary.txt')
        coverage_data = parse_coverage_file(coverage_file) if os.path.exists(coverage_file) else {}
        
        # Combine all non-BUSCO data
        other_data = {**quast_data, **merqury_qv_data, **merqury_completeness_data, **coverage_data}
        
        # Check both general and specific BUSCO directories
        busco_types = ['busco_general_pilon', 'busco_specific_pilon']
        
        for busco_type in busco_types:
            # Direct path to expected BUSCO file location
            busco_pattern = os.path.join(sample_path, 
                                        'best_assembly_info_and_QC', 
                                        busco_type, 
                                        'short_summary.specific.*.txt')
            
            busco_files = glob.glob(busco_pattern)
            
            if busco_files:
                busco_file = busco_files[0]  # Take first match
                best_source_file = os.path.join(sample_path, 'best_assembly_source.txt')
                
                filename = os.path.basename(busco_file)
                sample_id, lineage, reads_type, kmer_strategy, assembler = extract_sample_info(filename, best_source_file)
                
                if sample_id and lineage:
                    busco_results = parse_busco_file(busco_file)
                    if busco_results:
                        entry = {
                            'sample_id': sample_id,
                            'lineage': lineage,
                            'busco_type': busco_type.replace('busco_', '').replace('_pilon', ''),
                            'reads_type': reads_type,
                            'kmer_strategy': kmer_strategy,
                            'assembler': assembler,
                        }
                        # Add BUSCO results
                        entry.update(busco_results)
                        # Add all other metrics (QUAST, Merqury, Coverage)
                        entry.update(other_data)
                        
                        data.append(entry)
                        
                        if verbose:
                            print(f"✓ {sample_id} ({busco_type.replace('busco_', '').replace('_pilon', '')}) + All Metrics")
    
    return pd.DataFrame(data)


def create_wide_format_with_all_metrics(df):
    """
    Create wide format with one row per sample, including BUSCO, QUAST, Merqury, and Coverage metrics.
    """
    # Separate dataframes for general and specific
    df_general = df[df['busco_type'] == 'general'].copy()
    df_specific = df[df['busco_type'] == 'specific'].copy()
    
    # Define BUSCO metrics to include
    busco_metrics = ['complete_percent', 'single_copy_percent', 'duplicated_percent', 
                     'fragmented_percent', 'missing_percent', 'total_buscos',
                     'complete_count', 'single_copy_count', 'duplicated_count',
                     'fragmented_count', 'missing_count']
    
    # Use only the QUAST metrics that actually exist in your data
    quast_metrics = ['total_contigs>250bp', 'contigs>1000bp', 'contigs>5000bp', 'contigs>10000bp',
                     'contigs>25000bp', 'total_length>250bp', 'largest_contig', 'gc_percent', 'aUN',
                     'N50', 'N90', 'L50', 'L90', 'Ns_per_100kbp']
    
    # Define Merqury metrics
    merqury_metrics = ['error_percent', 'qv', 'percent_covered']
    
    # Define Coverage metrics
    coverage_metrics = ['mean_coverage', 'weighted_mean_coverage', 'median_coverage', 
                       'peak_coverage', 'min_coverage', 'std_coverage',
                       'mean_mapping_quality', 'median_mapping_quality',
                       'total_reads', 'mapped_reads', 'unmapped_reads',
                       'mapping_rate', 'properly_paired']
    
    # Add suffix to BUSCO metric columns
    for metric in busco_metrics:
        if metric in df_general.columns:
            df_general = df_general.rename(columns={metric: f'{metric}_general'})
        if metric in df_specific.columns:
            df_specific = df_specific.rename(columns={metric: f'{metric}_specific'})
    
    # Rename lineage columns
    df_general = df_general.rename(columns={'lineage': 'lineage_general'})
    df_specific = df_specific.rename(columns={'lineage': 'lineage_specific'})
    
    # Select columns to merge
    base_columns = ['sample_id', 'reads_type', 'kmer_strategy', 'assembler']  # Added kmer_strategy
    general_columns = (base_columns + ['lineage_general'] + 
                      [f'{m}_general' for m in busco_metrics if f'{m}_general' in df_general.columns])
    specific_columns = (['sample_id', 'lineage_specific'] + 
                       [f'{m}_specific' for m in busco_metrics if f'{m}_specific' in df_specific.columns] +
                      quast_metrics + merqury_metrics + coverage_metrics)  # All non-BUSCO metrics from specific df
    
    # Merge on sample_id
    result = df_general[general_columns].merge(
        df_specific[specific_columns], 
        on='sample_id', 
        how='outer'
    )
    
    return result




def main():
    """Main function to handle command line arguments and run the extraction."""
    parser = argparse.ArgumentParser(
        description='Extract comprehensive assembly metrics from genome assembly QC results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python extract_metrics.py /path/to/final_genomes_directory output_metrics.csv
  
  # Quiet mode (no progress messages)
  python extract_metrics.py /path/to/final_genomes_directory output_metrics.csv --quiet
        """
    )
    
    parser.add_argument('input_dir', 
                        help='Directory containing sample folders with assembly QC results')
    parser.add_argument('output_file', 
                        help='Output CSV file path')
    parser.add_argument('--quiet', '-q',
                        action='store_true',
                        help='Suppress progress messages')
    
    args = parser.parse_args()
    
    # Validate input directory
    if not os.path.exists(args.input_dir):
        print(f"Error: Input directory '{args.input_dir}' does not exist")
        sys.exit(1)
    
    if not os.path.isdir(args.input_dir):
        print(f"Error: '{args.input_dir}' is not a directory")
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(args.output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    try:
        # Collect all metrics
        if not args.quiet:
            print(f"Collecting metrics from: {args.input_dir}")
        
        df = collect_all_metrics(args.input_dir, verbose=not args.quiet)
        
        if df.empty:
            print("Warning: No data collected. Please check your input directory structure.")
            sys.exit(1)
        
        # Convert to wide format (one row per sample)
        if not args.quiet:
            print("Converting to wide format...")
        df_wide = create_wide_format_with_all_metrics(df)
        
        # Save results
        df_wide.to_csv(args.output_file, index=False)
        
        if not args.quiet:
            print(f"Results saved to: {args.output_file}")
            print(f"Final dataset shape: {df_wide.shape}")
            print(f"Samples processed: {df_wide['sample_id'].nunique() if 'sample_id' in df_wide.columns else 'N/A'}")
    
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()