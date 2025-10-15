#!/bin/bash

# This script is designed to save the output of the FSP assembly benchmarking pipeline easily.
# We will save: assembly fasta files, benchmarks, logs, QC results

# Function to show help
show_help() {
    cat << EOF
Usage: $0 -i <path_to_workflow_dir> -r <results_directory_name> -w <where_to_save> -n <name>

This script saves the output of the FSP assembly benchmarking pipeline.

Options:
  -i <path>   Path to the workflow root directory (FSP_assembly_benchmarking)
  -r <name>   Name of the directory containing the results of the workflow
  -w <path>   Path where the results should be saved
  -n <name>   Descriptive name for the run
  -h          Show this help message

Example:
  $0 -i /path/to/FSP_assembly_benchmarking -r results_trial2 -w /storage/final -n R1R2_default_kmer

EOF
}

# Initialise variables
workflow_path=""
results_directory_name=""
where_to_save=""
name=""

# getopts to set user-defined variables
while getopts "i:r:w:n:h" opt; do
  case $opt in
    i) workflow_path="$OPTARG"
       echo "The path to the workflow root directory: $OPTARG"
    ;;
    r) results_directory_name="$OPTARG"
       echo "The name of the directory containing the results: $OPTARG"
    ;;
    w) where_to_save="$OPTARG"
       echo "The path where the results should be saved: $OPTARG"
    ;;
    n) name="$OPTARG"
       echo "Chosen name for the run: $OPTARG"
    ;;
    h) show_help
       exit 0  # ← Exit after showing help
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
        show_help
        exit 1  # ← Exit with error for invalid options
    ;;
  esac
done

# save fasta files. All the fasta files are located in $results_directory_name/assemblies
mkdir -p "$where_to_save/$name"
cp -rL "$workflow_path/$results_directory_name/assemblies/" "$where_to_save/$name/"

# save the config file for the run
cp "$workflow_path/config/config.yml" "$where_to_save/$name/${name}_config.yml"

# save the log files for the run
cp -r "$workflow_path/logs/" "$where_to_save/$name/${name}_logs/"

# save the benchmark results for the run
cp -r "$workflow_path/benchmark/" "$where_to_save/$name/${name}_benchmark/"

# save the quast results for the run for each sample. The structure directory should be $workflow_path/$results_directory_name/$sample_name/quast

# Save QUAST results for each sample with sample name preserved
echo "Saving QUAST results..."
for sample_dir in "$workflow_path/$results_directory_name"/*/quast; do
    if [[ -d "$sample_dir" ]]; then
        # Extract sample name from path
        sample_name=$(basename "$(dirname "$sample_dir")")
        
        # Create descriptive directory name
        target_dir="$where_to_save/$name/QC/${sample_name}/quast"
        
        echo "  Copying $sample_name QUAST results to $target_dir"
        mkdir -p "$target_dir"
        cp -r "$sample_dir"/* "$target_dir/"
    fi
done


# Save BUSCO results - only short_summary.txt files
echo "Saving BUSCO summary files..."
for busco_dir in "$workflow_path/$results_directory_name"/*/busco_*; do
    if [[ -d "$busco_dir" ]]; then
        # Extract sample name and busco type
        sample_name=$(basename "$(dirname "$busco_dir")")
        busco_type=$(basename "$busco_dir")
        
        # Create descriptive directory name
        target_dir="$where_to_save/$name/QC/${sample_name}/${busco_type}"
        
        echo "  Copying $sample_name $busco_type summary files..."
        
        # Find and copy all short_summary.txt files (including assembler subdirectories)
        find "$busco_dir" -name "short_summary*.txt" | while read summary_file; do
            # Get the relative path from busco_dir to preserve structure
            relative_path="${summary_file#$busco_dir/}"
            target_file="$target_dir/$relative_path"
            
            # Create subdirectory if needed
            mkdir -p "$(dirname "$target_file")"
            
            # Copy the summary file
            cp "$summary_file" "$target_file"
            echo "    Copied: $relative_path"
        done
    fi
done


# Save MerquryFK results - key files only
echo "Saving MerquryFK results..."
for merqury_dir in "$workflow_path/$results_directory_name"/*/merquryfk; do
    if [[ -d "$merqury_dir" ]]; then
        sample_name=$(basename "$(dirname "$merqury_dir")")
        
        echo "  Copying $sample_name MerquryFK results..."
        
        # Iterate through each assembler subdirectory
        for assembler_dir in "$merqury_dir"/*/; do
            if [[ -d "$assembler_dir" ]]; then
                assembler_name=$(basename "$assembler_dir")
                target_dir="$where_to_save/$name/QC/${sample_name}/merquryfk/${assembler_name}"
                
                mkdir -p "$target_dir"
                
                # Copy key MerquryFK files (avoid large intermediate files)
                echo "    Copying $assembler_name results..."
                
                # Copy QV files (quality value - most important)
                find "$assembler_dir" -name "*.qv" -exec cp {} "$target_dir/" \;
                
                # Copy completeness stats
                find "$assembler_dir" -name "merquryfk.completeness.stats" -exec cp {} "$target_dir/" \;
                
                # Copy plots (PNG files)
                find "$assembler_dir" -name "*.png" -exec cp {} "$target_dir/" \;
                
                # Copy BED files (optional - remove if too large)
                find "$assembler_dir" -name "*.bed" -exec cp {} "$target_dir/" \;
                
                echo "      Copied $assembler_name MerquryFK files"
            fi
        done
    fi
done


echo "Results saved successfully to: $where_to_save/$name"