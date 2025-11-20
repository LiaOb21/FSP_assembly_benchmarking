#!/bin/bash

# This script is designed to save the output of the FSP assembly benchmarking pipeline easily.
# We will save: assembly fasta files, benchmarks, logs, QC results

# Function to show help
show_help() {
    cat << EOF
Usage: $0 -i <path_to_workflow_dir> -r <results_directory_name> -w <where_to_save> -b <batch_number> -t <reads_type> -k <kmer_strategy> -n <run_id>

This script saves the output of the FSP assembly benchmarking pipeline.

Options:
  -i <path>   Path to the workflow root directory (FSP_assembly_benchmarking)
  -r <name>   Name of the directory containing the results of the workflow
  -w <path>   Path where the results should be saved
  -b <batch>  Batch number (e.g., EG1, batch2, B01, etc.). If you used multiple batches, use all the identifiers (e.g., EG1_EG2). If you don't have batches, use whatever you like in this field.
  -t <type>   Reads used: R1R2 or merged
  -k <kmer>   Kmer strategy: manual, kmergenie, or reads_length
  -n <run_id> Run identifier (e.g., run1, subset_A, test, trial1, etc.). This should be something that describes your run, for example number of samples + kmer distribution (e.g., 13samples_SP, 34samples_DD, etc. SP = single peak, DD = double diploid peak). Allows to distinguish between multiple runs with the same batch, reads type, and kmer strategy.
  -h          Show this help message

The run name will be automatically generated as: <batch>_<reads_type>_<kmer_strategy>_<run_id>

Examples:
  $0 -i /path/to/FSP_assembly_benchmarking -r results -w /storage/final -b batch1 -t R1R2 -k manual -n run1
  → Creates: batch1_R1R2_manual_run1

  $0 -i . -r results -w ./saved_output -b B02 -t merged -k kmergenie -n subset_A
  → Creates: B02_merged_kmergenie_subset_A

EOF
}

# Function to validate expected file structures
validate_sample_structure() {
    local workflow_path="$1"
    local results_dir="$2"
    local validation_errors=0

    echo "Validating sample directory structures..."

    # Get sample list (excluding fqreads)
    local sample_count=0
    for sample_path in "$workflow_path/$results_dir"/*/; do
        if [[ -d "$sample_path" && "$(basename "$sample_path")" != "fqreads" ]]; then
            local sample_name=$(basename "$sample_path")
            ((sample_count++))
            
            echo "  Validating sample directory structure: $sample_name"
            
            # Check required directories
            local required_dirs=(
                "assemblies"
                "quast" 
                "busco_general"
                "busco_specific"
                "merquryfk"
                "best_assembly"
                "best_assembly_qc"
            )
            
            for dir in "${required_dirs[@]}"; do
                if [[ ! -d "$sample_path/$dir" ]]; then
                    echo "    🍄 Missing required directory: $dir"
                    ((validation_errors++))
                else
                    echo "    ✓ Found directory: $dir"
                fi
            done
            
            # Check for assembly files
            if [[ -d "$sample_path/assemblies" ]]; then
                if ls "$sample_path/assemblies"/*.fa 1> /dev/null 2>&1; then
                    local fa_count=$(ls "$sample_path/assemblies"/*.fa | wc -l)
                    echo "    ✓ Found $fa_count assembly files (including best assembly improved with pypolca)"
                else
                    echo "    🍄 No .fa files found in assemblies directory"
                    ((validation_errors++))
                fi
            fi
            
            # Check for QUAST reports
            if [[ -d "$sample_path/quast" ]]; then
                if ls "$sample_path/quast"/*report* 1> /dev/null 2>&1; then
                    echo "    ✓ Found QUAST reports"
                else
                    echo "    🍄 Warning: No QUAST reports found"
                    ((validation_errors++))
                fi
            fi
            
            # Check for BUSCO results
            for busco_type in "busco_general" "busco_specific"; do
                if [[ -d "$sample_path/$busco_type" ]]; then
                    if find "$sample_path/$busco_type" -name "short_summary.*.txt" | grep -q .; then
                        echo "    ✓ Found $busco_type results"
                    else
                        echo "    🍄 Warning: No BUSCO summaries found in $busco_type"
                        ((validation_errors++))
                    fi
                fi
            done

            # Check for merquryfk results
            if [[ -d "$sample_path/merquryfk" ]]; then
                if ls "$sample_path/merquryfk"/*/*{qv,stats,bed,png} 1> /dev/null 2>&1; then
                    echo "    ✓ Found MerquryFK results"
                else
                    echo "    🍄 Warning: No MerquryFK results found"
                    ((validation_errors++))
                fi
            fi

            # Check for best assembly info
            if [[ -d "$sample_path/best_assembly" ]]; then
                if ls "$sample_path/best_assembly"/*.txt 1> /dev/null 2>&1; then
                    echo "    ✓ Found best assembly info"
                else
                    echo "    🍄 Warning: No best assembly info found"
                    ((validation_errors++))
                fi
            fi            
            
            # Check for pypolca vcf file
            if [[ -d "$sample_path/best_assembly" ]]; then
                if ls "$sample_path/best_assembly"/pypolca/*.vcf 1> /dev/null 2>&1; then
                    echo "    ✓ Found pypolca vcf file"
                else
                    echo "    🍄 Warning: No pypolca vcf file found"
                    ((validation_errors++))
                fi
            fi  

            # Check for best assembly QUAST reports
            if [[ -d "$sample_path/best_assembly_qc/quast_pypolca" ]]; then
                if ls "$sample_path/best_assembly_qc/quast_pypolca"/*report* 1> /dev/null 2>&1; then
                    echo "    ✓ Found best assembly QUAST reports"
                else
                    echo "    🍄 Warning: No best assembly QUAST reports found"
                    ((validation_errors++))
                fi
            else
                echo "    🍄 Warning: Best assembly QUAST directory not found"
                ((validation_errors++))
            fi

            # Check for best assembly BUSCO results
            for busco_type in "busco_general_pypolca" "busco_specific_pypolca"; do
                if [[ -d "$sample_path/best_assembly_qc/$busco_type" ]]; then
                    if find "$sample_path/best_assembly_qc/$busco_type" -name "short_summary.*.txt" | grep -q .; then
                        echo "    ✓ Found best assembly $busco_type results"
                    else
                        echo "    🍄 Warning: No best assembly BUSCO summaries found in $busco_type"
                        ((validation_errors++))
                    fi
                else
                    echo "    🍄 Warning: Best assembly BUSCO directory $busco_type not found"
                    ((validation_errors++))
                fi
            done
            
            # Check for best assembly merquryfk results
            if [[ -d "$sample_path/best_assembly_qc/merquryfk_pypolca" ]]; then
                if ls "$sample_path/best_assembly_qc/merquryfk_pypolca"/*{qv,stats,bed,png} 1> /dev/null 2>&1; then
                    echo "    ✓ Found best assembly MerquryFK results"
                else
                    echo "    🍄 Warning: No best assembly MerquryFK results found"
                    ((validation_errors++))
                fi
            else
                echo "    🍄 Warning: Best assembly MerquryFK directory not found"
                ((validation_errors++))
            fi

            # Check for best assembly alignments results
            if [[ -d "$sample_path/best_assembly_qc/samtools_pypolca" ]]; then
                if ls "$sample_path/best_assembly_qc/samtools_pypolca"/*{sorted.bam,.txt} 1> /dev/null 2>&1; then
                    echo "    ✓ Found best assembly alignments results"
                else
                    echo "    🍄 Warning: No best assembly alignments results found"
                    ((validation_errors++))
                fi
            else
                echo "    🍄 Warning: Best assembly samtools_pypolca directory not found"
                ((validation_errors++))
            fi
            # Check for best assembly coverage_viz results
            if [[ -d "$sample_path/best_assembly_qc/coverage_viz_pypolca" ]]; then
                if ls "$sample_path/best_assembly_qc/coverage_viz_pypolca"/*{.png,.txt} 1> /dev/null 2>&1; then
                    echo "    ✓ Found best assembly coverage_viz results"
                else
                    echo "    🍄 Warning: No best assembly coverage_viz results found"
                    ((validation_errors++))
                fi
            else
                echo "    🍄 Warning: Best assembly coverage_viz_pypolca directory not found"
                ((validation_errors++))
            fi
        fi
    done
    
    # Final validation summary
    if [[ $sample_count -eq 0 ]]; then
        echo "🍄 ERROR: No sample directories found!" >&2
        return 1
    fi
    
    if [[ $validation_errors -gt 0 ]]; then
        echo "🍄 ERROR: Found $validation_errors validation errors across samples!" >&2
        echo "Some expected files/directories are missing. Continue anyway? (y/N)"
        read -r response
        if [[ ! "$response" =~ ^[Yy]$ ]]; then
            echo "Operation cancelled."
            return 1
        fi
    else
        echo "✓ All validation checks passed for $sample_count samples!"
        echo "--------------------------------------------------------------"
    fi
    
    return 0
}

# Initialise variables
workflow_path=""
results_directory_name=""
where_to_save=""
batch_number=""
reads_type=""
kmer_strategy=""
run_id=""

# getopts to set user-defined variables
while getopts "i:r:w:b:t:k:n:h" opt; do
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
    b) batch_number="$OPTARG"
       echo "Batch number: $OPTARG"
    ;;
    t) reads_type="$OPTARG"
       # Validate reads type
       if [[ "$OPTARG" != "R1R2" && "$OPTARG" != "merged" ]]; then
           echo "Error: Reads type must be 'R1R2' or 'merged'" >&2
           exit 1
       fi
       echo "Reads type: $OPTARG"
    ;;
    k) kmer_strategy="$OPTARG"
       # Validate kmer strategy
       if [[ "$OPTARG" != "manual" && "$OPTARG" != "kmergenie" && "$OPTARG" != "reads_length" ]]; then
           echo "Error: Kmer strategy must be 'manual', 'kmergenie', or 'reads_length'" >&2
           exit 1
       fi
       echo "Kmer strategy: $OPTARG"
    ;;
    n) run_id="$OPTARG"
       echo "Run identifier: $OPTARG"
    ;;
    h) show_help
       exit 0
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
        show_help
        exit 1
    ;;
  esac
done

# Check that all required parameters are provided
if [[ -z "$workflow_path" || -z "$results_directory_name" || -z "$where_to_save" || -z "$batch_number" || -z "$reads_type" || -z "$kmer_strategy" || -z "$run_id" ]]; then
    echo "Error: All parameters are required!" >&2
    echo ""
    show_help
    exit 1
fi

# Check if the workflow path exists
if [[ ! -d "$workflow_path" ]]; then
    echo "🍄 Error: Workflow path '$workflow_path' does not exist!" >&2
    exit 1
fi

# Check if the results directory exists
if [[ ! -d "$workflow_path/$results_directory_name" ]]; then
    echo "🍄 Error: Results directory '$workflow_path/$results_directory_name' does not exist!" >&2
    exit 1
fi

#Check if the benchmark directory exists if not, warn the user and ask for confirmation to continue
if [[ ! -d "$workflow_path/benchmark" ]]; then
    echo "🍄 Warning: Benchmark directory '$workflow_path/benchmark' does not exist!" >&2
    echo "Continue without saving benchmark results? (y/N)"
    read -r response
    if [[ ! "$response" =~ ^[Yy]$ ]]; then
        echo "Operation cancelled."
        exit 1
    fi
fi

#Check if the logs directory exists if not, warn the user and ask for confirmation to continue
if [[ ! -d "$workflow_path/logs" ]]; then
    echo "🍄 Warning: Logs directory '$workflow_path/logs' does not exist!" >&2
    echo "Continue without saving log files? (y/N)"
    read -r response
    if [[ ! "$response" =~ ^[Yy]$ ]]; then
        echo "Operation cancelled."
        exit 1
    fi
fi

# Check if the config file exists if not, warn the user and ask for confirmation to continue
if [[ ! -f "$workflow_path/config/config.yml" ]]; then
    echo "🍄 Warning: Config file '$workflow_path/config/config.yml' does not exist!" >&2
    echo "Continue without saving config file? (y/N)"
    read -r response
    if [[ ! "$response" =~ ^[Yy]$ ]]; then
        echo "Operation cancelled."
        exit 1
    fi
fi

# Generate the run name automatically
name="${batch_number}_${reads_type}_${kmer_strategy}_${run_id}"
echo ""
echo "Generated run name: $name"
echo ""

# Check if a directory with the same name as the target directory already exists to avoid accidental overwriting
if [[ -d "$where_to_save/$name" ]]; then
    echo "🍄 WARNING: Directory '$where_to_save/$name' already exists!"
    echo "This will overwrite existing files. Continue? (y/N)"
    read -r response
    if [[ ! "$response" =~ ^[Yy]$ ]]; then
        echo "Operation cancelled."
        exit 1
    fi
fi

# Validate sample structures before starting
if ! validate_sample_structure "$workflow_path" "$results_directory_name"; then
    exit 1
fi

# Create the target directory
echo "Creating target directory..."
mkdir -p "$where_to_save/$name"

# Produce a simple text file containing the list of the samples used in the run
echo "Saving sample list..."
ls "$workflow_path/$results_directory_name" | grep -v "fqreads" > "$where_to_save/$name/${name}_sample_list.txt"

# save fasta files. All the fasta files produced by the workflow are located in $results_directory_name/$sample_name/assemblies (including the best assembly improved with pilon)
echo "Saving fasta files..."
for sample_dir in "$workflow_path/$results_directory_name"/*/assemblies; do
    if [[ -d "$sample_dir" ]]; then
        # Extract sample name from path (e.g., results/048ds/assemblies -> 048ds)
        sample_name=$(basename "$(dirname "$sample_dir")")
        
        # Create target directory maintaining sample structure
        target_dir="$where_to_save/$name/$sample_name/assemblies"
        
        echo "Copying $sample_name assemblies to $target_dir"
        mkdir -p "$target_dir"
        cp -rL "$sample_dir"/*.fa "$target_dir/"
    fi
done

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
        target_dir="$where_to_save/$name/${sample_name}/QC_all_drafts/quast"
        
        echo "  Copying $sample_name QUAST results to $target_dir"
        mkdir -p "$target_dir"
        cp -r "$sample_dir"/*report* "$target_dir/"
    fi
done


# Save BUSCO results - only short_summary.txt files and full_table.tsv for downstream analysis
echo "Saving BUSCO summary files..."
for busco_dir in "$workflow_path/$results_directory_name"/*/busco_{general,specific}; do
    if [[ -d "$busco_dir" ]]; then
        # Extract sample name and busco type
        sample_name=$(basename "$(dirname "$busco_dir")")
        busco_type=$(basename "$busco_dir")
        
        # Create descriptive directory name
        target_dir="$where_to_save/$name/${sample_name}/QC_all_drafts/${busco_type}"
        
        echo "  Copying $sample_name $busco_type summary files..."
        
        # Process each assembler subdirectory individually to avoid duplicates
        for assembler_subdir in "$busco_dir"/*/*/; do
            if [[ -d "$assembler_subdir" ]]; then
                assembler_name=$(basename "$assembler_subdir" .fa)
                
                # Create target directory
                target_dir="$where_to_save/$name/${sample_name}/QC_all_drafts/${busco_type}/${assembler_name}"
                mkdir -p "$target_dir"
                
                # Copy only the main short_summary file (not from subdirectories)
                short_summary_file="$assembler_subdir/short_summary.specific."*".txt"
                if ls $short_summary_file 1> /dev/null 2>&1; then
                    cp $short_summary_file "$target_dir/"
                    echo "    Copied $assembler_name BUSCO summary"
                fi
                
                # Also copy full_table.tsv for downstream analysis
                full_table_file="$assembler_subdir/run_*/full_table.tsv"
                if ls $full_table_file 1> /dev/null 2>&1; then
                    cp $full_table_file "$target_dir/${assembler_name}_full_table.tsv"
                    echo "    Copied $assembler_name BUSCO full table"
                fi
            fi
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
                target_dir="$where_to_save/$name/${sample_name}/QC_all_drafts/merquryfk/${assembler_name}"
                
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


# save best assembly information and QC
echo "Saving best assembly info..."
for sample_dir in "$workflow_path/$results_directory_name"/*/best_assembly; do
    if [[ -d "$sample_dir" ]]; then
        # Extract sample name from path
        sample_name=$(basename "$(dirname "$sample_dir")")
        
        # Create descriptive directory name
        target_dir="$where_to_save/$name/${sample_name}/best_assembly_info_and_QC/info"

        echo "  Copying $sample_name best assembly info to $target_dir"
        mkdir -p "$target_dir"
        cp -r "$sample_dir"/*txt "$target_dir/"
        cp -r "$sample_dir"/pypolca/*vcf "$target_dir/"
        cp -r "$sample_dir"/pypolca/*report "$target_dir/"
    fi
done


# Save QUAST results for each sample with sample name preserved
echo "Saving QUAST results for best assemblies improved with pypolca..."
for sample_dir in "$workflow_path/$results_directory_name"/*/best_assembly_qc/quast_pypolca; do
    if [[ -d "$sample_dir" ]]; then
        # Extract sample name from path
        sample_name=$(basename "$(dirname "$(dirname "$sample_dir")")")
        
        # Create descriptive directory name
        target_dir="$where_to_save/$name/${sample_name}/best_assembly_info_and_QC/quast_pypolca"
        
        echo "  Copying $sample_name QUAST results to $target_dir"
        mkdir -p "$target_dir"
        cp -r "$sample_dir"/*report* "$target_dir/"
    fi
done


# Save BUSCO results - only short_summary.txt files and full_table.tsv for downstream analysis
echo "Saving best assembly BUSCO summary files..."
for busco_dir in "$workflow_path/$results_directory_name"/*/best_assembly_qc/busco_*; do
    if [[ -d "$busco_dir" ]]; then
        # Extract sample name and busco type
        sample_name=$(basename "$(dirname "$(dirname "$busco_dir")")")
        busco_type=$(basename "$busco_dir")
        
        # Create descriptive directory name
        target_dir="$where_to_save/$name/${sample_name}/best_assembly_info_and_QC/${busco_type}"
        mkdir -p "$target_dir"
        
        echo "  Copying $sample_name $busco_type summary files..."
        
        # Copy the main short_summary file (directly in busco_dir, no assembler subdirs)
        short_summary_file="$busco_dir/*/short_summary.specific."*".txt"
        if ls $short_summary_file 1> /dev/null 2>&1; then
            cp $short_summary_file "$target_dir/"
            echo "    Copied best assembly BUSCO summary"
        fi

        # Also copy full_table.tsv for downstream analysis
        full_table_file="$busco_dir/*/run_*/full_table.tsv"
        if ls $full_table_file 1> /dev/null 2>&1; then
            cp $full_table_file "$target_dir/${sample_name}_best_assembly_full_table.tsv"
            echo "    Copied $busco_type BUSCO full table"
        fi
    fi
done


# Save MerquryFK results - key files only
echo "Saving best assembly MerquryFK results..."
for merqury_dir in "$workflow_path/$results_directory_name"/*/best_assembly_qc/merquryfk_pypolca; do
    if [[ -d "$merqury_dir" ]]; then
        # Extract sample name  
        sample_name=$(basename "$(dirname "$(dirname "$merqury_dir")")")

        # Create target directory
        target_dir="$where_to_save/$name/${sample_name}/best_assembly_info_and_QC/merquryfk_pypolca"
        mkdir -p "$target_dir"
        
        echo "  Copying $sample_name best assembly MerquryFK results..."
        
        # Copy key MerquryFK files (files are directly in merqury_dir, no assembler subdirs)
        
        # Copy QV files (quality value - most important)
        find "$merqury_dir" -maxdepth 1 -name "*.qv" -exec cp {} "$target_dir/" \;
        
        # Copy completeness stats
        find "$merqury_dir" -maxdepth 1 -name "merquryfk.completeness.stats" -exec cp {} "$target_dir/" \;
        
        # Copy plots (PNG files)
        find "$merqury_dir" -maxdepth 1 -name "*.png" -exec cp {} "$target_dir/" \;
        
        # Copy BED files (optional - remove if too large)
        find "$merqury_dir" -maxdepth 1 -name "*.bed" -exec cp {} "$target_dir/" \;
        
        echo "    Copied best assembly MerquryFK files"
    fi
done

# Save alignments of reads to best assemblies improved with pypolca (BAM files and indexes)
echo "Saving alignments and stats results for best assemblies improved with pypolca..."
for sample_dir in "$workflow_path/$results_directory_name"/*/best_assembly_qc/samtools_pypolca; do
    if [[ -d "$sample_dir" ]]; then
        # Extract sample name from path
        sample_name=$(basename "$(dirname "$(dirname "$sample_dir")")")
        
        # Create descriptive directory name
        target_dir="$where_to_save/$name/${sample_name}/best_assembly_info_and_QC/samtools_pypolca"

        echo "  Copying $sample_name samtools_pypolca results to $target_dir"
        mkdir -p "$target_dir"
        cp -r "$sample_dir"/*sorted.bam* "$target_dir/"
        cp -r "$sample_dir"/*.txt "$target_dir/"
    fi
done

# Save coverage visualisation results for best assemblies improved with pypolca
echo "Saving coverage visualisation results for best assemblies improved with pypolca..."
for sample_dir in "$workflow_path/$results_directory_name"/*/best_assembly_qc/coverage_viz_pypolca; do
    if [[ -d "$sample_dir" ]]; then
        # Extract sample name from path
        sample_name=$(basename "$(dirname "$(dirname "$sample_dir")")")
        
        # Create descriptive directory name
        target_dir="$where_to_save/$name/${sample_name}/best_assembly_info_and_QC/coverage_viz_pypolca"

        echo "  Copying $sample_name coverage_viz_pypolca results to $target_dir"
        mkdir -p "$target_dir"
        cp -r "$sample_dir"/*png "$target_dir/"
        cp -r "$sample_dir"/*.txt "$target_dir/"
    fi
done


echo "Results saved successfully to: $where_to_save/$name"