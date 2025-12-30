#!/bin/bash

# This script is designed to save the output of the FSP assembly benchmarking pipeline easily.
# We will save: assembly fasta files, benchmarks, logs, QC results

# Function to show help
show_help() {
    cat << EOF
Usage: $0 -i <path_to_workflow_dir> -r <results_directory_name> -w <where_to_save> -b <batch_number> -n <run_id> [-f <final_best_dir>]

This script saves the output of the FSP assembly benchmarking pipeline.

Options:
  -i <path>   Path to the workflow root directory (FSP_assembly_benchmarking)
  -r <name>   Name of the directory containing the results of the workflow
  -w <path>   Path where the results should be saved
  -b <batch>  Batch number (e.g., EG1, batch2, B01, etc.). If you used multiple batches, use all the identifiers (e.g., EG1_EG2). If you don't have batches, use whatever you like in this field.
  -n <run_id> Run identifier (e.g., run1, subset_A, test, trial1, etc.). This should be something that describes your run.
  -f <path>   [Optional] Full path where to save final best assemblies. If not provided, they will be saved in <where_to_save>/<run_name>/FINAL_BEST_ASSEMBLIES
  -h          Show this help message

The run name will be automatically generated as: <batch>_<run_id>

Examples:
  $0 -i /path/to/FSP_assembly_benchmarking -r results -w /storage/final -b batch1 -n run1
  → Creates: /storage/final/batch1_run1/
  → Final best assemblies: /storage/final/batch1_run1/FINAL_BEST_ASSEMBLIES/

  $0 -i . -r results_3 -w ./saved_output -b EG1 -n 48samples -f /storage/FINAL_ASSEMBLIES/EG1_48samples
  → Creates: ./saved_output/EG1_48samples/
  → Final best assemblies: /storage/FINAL_ASSEMBLIES/EG1_48samples/

EOF
}


# Function to validate expected file structures
validate_sample_structure() {
    local workflow_path="$1"
    local results_dir="$2"
    local validation_errors=0

    echo "Validating results directory structure..."

    # Get sample list from assemblies directory
    local sample_count=0
    for sample_path in "$workflow_path/$results_dir/assemblies"/*/; do
        if [[ -d "$sample_path" ]]; then
            local sample_name=$(basename "$sample_path")
            ((sample_count++))
            
            echo "  Validating files for sample: $sample_name"
            
            # Check for assembly files
            if ls "$sample_path"/*.fa 1> /dev/null 2>&1; then
                local fa_count=$(ls "$sample_path"/*.fa | wc -l)
                echo "    ✓ Found $fa_count assembly files"
            else
                echo "    🍄 No .fa files found for $sample_name"
                ((validation_errors++))
            fi
            
            # Check for QUAST reports
            if [[ -d "$workflow_path/$results_dir/quast/$sample_name" ]]; then
                if ls "$workflow_path/$results_dir/quast/$sample_name"/*report* 1> /dev/null 2>&1; then
                    echo "    ✓ Found QUAST reports"
                else
                    echo "    🍄 Warning: No QUAST reports found"
                    ((validation_errors++))
                fi
            else
                echo "    🍄 Warning: QUAST directory not found"
                ((validation_errors++))
            fi
            
            # Check for BUSCO results across reads_types and strategies
            local busco_found=0
            for reads_type in R1R2 merged; do
                for strategy in manual kmergenie reads_length; do
                    if [[ -d "$workflow_path/$results_dir/$reads_type/$strategy/$sample_name" ]]; then
                        for busco_type in busco_general busco_specific; do
                            if [[ -d "$workflow_path/$results_dir/$reads_type/$strategy/$sample_name/$busco_type" ]]; then
                                ((busco_found++))
                            fi
                        done
                    fi
                done
            done
            
            if [[ $busco_found -gt 0 ]]; then
                echo "    ✓ Found BUSCO results in $busco_found locations"
            else
                echo "    🍄 Warning: No BUSCO results found"
                ((validation_errors++))
            fi

            # Check for best assembly info
            if [[ -d "$workflow_path/$results_dir/best_assembly/$sample_name" ]]; then
                if ls "$workflow_path/$results_dir/best_assembly/$sample_name"/*.txt 1> /dev/null 2>&1; then
                    echo "    ✓ Found best assembly info"
                else
                    echo "    🍄 Warning: No best assembly info found"
                    ((validation_errors++))
                fi
            else
                echo "    🍄 Warning: Best assembly directory not found"
                ((validation_errors++))
            fi
            
            # Check for best assembly QC
            if [[ -d "$workflow_path/$results_dir/best_assembly_qc/$sample_name" ]]; then
                echo "    ✓ Found best assembly QC directory"
            else
                echo "    🍄 Warning: Best assembly QC directory not found"
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
        echo "🍄 Found $validation_errors validation warnings across samples!"
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
run_id=""
final_best_dir="" 

# getopts to set user-defined variables
while getopts "i:r:w:b:n:f:h" opt; do
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
    n) run_id="$OPTARG"
       echo "Run identifier: $OPTARG"
    ;;
    f) final_best_dir="$OPTARG"
       echo "Final best assemblies will be saved to: $OPTARG"
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
if [[ -z "$workflow_path" || -z "$results_directory_name" || -z "$where_to_save" || -z "$batch_number" || -z "$run_id" ]]; then
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

# Generate the run name automatically (simplified)
name="${batch_number}_${run_id}"
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

# Create the target directory
echo "Creating target directory..."
mkdir -p "$where_to_save/$name"

# Set up logging - redirect all output to both terminal and log file
log_file="$where_to_save/$name/${name}_save_output.log"
exec > >(tee -a "$log_file") 2>&1

echo "=============================================================="
echo "Logging to: $log_file"
echo "Script started at: $(date)"
echo "=============================================================="
echo ""

# Validate sample structures before starting
if ! validate_sample_structure "$workflow_path" "$results_directory_name"; then
    exit 1
fi

# Produce a simple text file containing the list of the samples used in the run
echo "Saving sample list..."
# Get sample names from assemblies directory
for sample_path in "$workflow_path/$results_directory_name/assemblies"/*/; do
    if [[ -d "$sample_path" ]]; then
        basename "$sample_path"
    fi
done > "$where_to_save/$name/${name}_sample_list.txt"

# save fasta files. All the fasta files are located in $workflow_path/$results_directory_name/assemblies/{sample}/
echo "Saving fasta files..."
for sample_dir in "$workflow_path/$results_directory_name/assemblies"/*; do
    if [[ -d "$sample_dir" ]]; then
        # Extract sample name from path (e.g., assemblies/048ds -> 048ds)
        sample_name=$(basename "$sample_dir")
        
        # Create target directory maintaining sample structure
        target_dir="$where_to_save/$name/$sample_name/assemblies"
        
        echo "  Copying $sample_name assemblies to $target_dir"
        mkdir -p "$target_dir"
        cp -L "$sample_dir"/*.fa "$target_dir/" 2>/dev/null || echo "    Warning: No .fa files found for $sample_name"
    fi
done

# save the config file for the run
echo "Saving config file..."
cp "$workflow_path/config/config.yml" "$where_to_save/$name/${name}_config.yml"

# save the log files for the run
echo "Saving log files..."
cp -r "$workflow_path/logs/" "$where_to_save/$name/${name}_logs/"

# save the benchmark results for the run
echo "Saving benchmark results..."
cp -r "$workflow_path/benchmark/" "$where_to_save/$name/${name}_benchmark/"

# save the quast results for the run for each sample. The structure directory should be $workflow_path/$results_directory_name/$sample_name/quast

# Save QUAST results for each sample with sample name preserved
echo "Saving QUAST results..."
for sample_dir in "$workflow_path/$results_directory_name/quast"/*; do
    if [[ -d "$sample_dir" ]]; then
        # Extract sample name from path (e.g., quast/048ds -> 048ds)
        sample_name=$(basename "$sample_dir")
        
        # Create descriptive directory name
        target_dir="$where_to_save/$name/${sample_name}/QC_all_drafts/quast"
        
        echo "  Copying $sample_name QUAST results to $target_dir"
        mkdir -p "$target_dir"
        
        # Copy report files
        if ls "$sample_dir"/*report.tsv 1> /dev/null 2>&1; then
            cp "$sample_dir"/*report.tsv "$target_dir/" 2>/dev/null
        else
            echo "    Warning: No QUAST reports found for $sample_name"
        fi
    fi
done


# Save BUSCO results - only short_summary.txt files and full_table.tsv for downstream analysis
echo "Saving BUSCO summary files..."

# Iterate through reads_types and strategies
for reads_type in R1R2 merged; do
    for strategy in manual kmergenie reads_length; do
        # Check if this combination exists
        strategy_dir="$workflow_path/$results_directory_name/$reads_type/$strategy"
        if [[ ! -d "$strategy_dir" ]]; then
            continue
        fi
        
        # Iterate through samples
        for sample_dir in "$strategy_dir"/*; do
            if [[ ! -d "$sample_dir" ]]; then
                continue
            fi
            
            sample_name=$(basename "$sample_dir")
            
            # Process both busco_general and busco_specific
            for busco_type in busco_general busco_specific; do
                busco_type_dir="$sample_dir/$busco_type"
                
                if [[ ! -d "$busco_type_dir" ]]; then
                    continue
                fi
                
                echo "  Processing $sample_name $busco_type ($reads_type/$strategy)..."
                
                # Iterate through assembler subdirectories
                for assembler_parent_dir in "$busco_type_dir"/*; do
                    if [[ ! -d "$assembler_parent_dir" ]]; then
                        continue
                    fi
                    
                    assembler_name=$(basename "$assembler_parent_dir")
                    
                    # BUSCO creates a nested BUSCO_* subdirectory
                    # Find the actual BUSCO output directory
                    busco_dir=$(find "$assembler_parent_dir" -mindepth 1 -maxdepth 1 -type d -name "BUSCO_*" | head -1)
                    
                    if [[ -z "$busco_dir" || ! -d "$busco_dir" ]]; then
                        echo "    Warning: No BUSCO output directory found in $assembler_parent_dir"
                        continue
                    fi
                    
                    # Create unique target directory including reads_type and strategy
                    target_dir="$where_to_save/$name/${sample_name}/QC_all_drafts/${busco_type}/${reads_type}_${strategy}_${assembler_name}"
                    mkdir -p "$target_dir"
                    
                    # Copy short_summary file - use find
                    if find "$busco_dir" -maxdepth 1 -name "short_summary.specific.*.txt" -type f -print -quit | grep -q .; then
                        find "$busco_dir" -maxdepth 1 -name "short_summary.specific.*.txt" -type f -exec cp {} "$target_dir/" \;
                        echo "    Copied ${reads_type}_${strategy}_${assembler_name} BUSCO summary"
                    fi
                    
                    # Copy full_table.tsv - use find
                    if find "$busco_dir" -path "*/run_*/full_table.tsv" -type f -print -quit | grep -q .; then
                        find "$busco_dir" -path "*/run_*/full_table.tsv" -type f -exec cp {} "$target_dir/${assembler_name}_full_table.tsv" \;
                        echo "    Copied ${reads_type}_${strategy}_${assembler_name} BUSCO full table"
                    fi
                done
            done
        done
    done
done

# Save MerquryFK results - key files only
echo "Saving MerquryFK results..."

# Iterate through reads_types and strategies
for reads_type in R1R2 merged; do
    for strategy in manual kmergenie reads_length; do
        # Check if this combination exists
        strategy_dir="$workflow_path/$results_directory_name/$reads_type/$strategy"
        if [[ ! -d "$strategy_dir" ]]; then
            continue
        fi
        
        # Iterate through samples
        for sample_dir in "$strategy_dir"/*; do
            if [[ ! -d "$sample_dir" ]]; then
                continue
            fi
            
            sample_name=$(basename "$sample_dir")
            
            # Check if merquryfk directory exists
            merqury_dir="$sample_dir/merquryfk"
            if [[ ! -d "$merqury_dir" ]]; then
                continue
            fi
            
            echo "  Processing $sample_name MerquryFK results ($reads_type/$strategy)..."
            
            # Iterate through assembler subdirectories
            for assembler_dir in "$merqury_dir"/*; do
                if [[ ! -d "$assembler_dir" ]]; then
                    continue
                fi
                
                assembler_name=$(basename "$assembler_dir")
                
                # Create unique target directory including reads_type and strategy
                target_dir="$where_to_save/$name/${sample_name}/QC_all_drafts/merquryfk/${reads_type}_${strategy}_${assembler_name}"
                mkdir -p "$target_dir"
                
                echo "    Copying ${reads_type}_${strategy}_${assembler_name} MerquryFK results..."
                
                # Copy QV files (quality value - most important)
                find "$assembler_dir" -name "*.qv" -exec cp {} "$target_dir/" \;
                
                # Copy completeness stats
                find "$assembler_dir" -name "merquryfk.completeness.stats" -exec cp {} "$target_dir/" \;
                
                # Copy plots (PNG files)
                find "$assembler_dir" -name "*ln.png" -exec cp {} "$target_dir/" \;
                
                # Copy BED files (optional - remove if too large)
                find "$assembler_dir" -name "*.bed" -exec cp {} "$target_dir/" \;
                
                echo "      Copied ${reads_type}_${strategy}_${assembler_name} MerquryFK files"
            done
        done
    done
done

# save best assembly information and QC
echo "Saving best assembly info..."
for sample_dir in "$workflow_path/$results_directory_name/best_assembly"/*; do
    if [[ -d "$sample_dir" ]]; then
        # Extract sample name from path (best_assembly/048ds -> 048ds)
        sample_name=$(basename "$sample_dir")
        
        # Create descriptive directory name
        target_dir="$where_to_save/$name/${sample_name}/best_assembly_info_and_QC/info"

        echo "  Copying $sample_name best assembly info to $target_dir"
        mkdir -p "$target_dir"
        
        # Copy txt files (best_assembly.txt, etc.)
        if ls "$sample_dir"/*.txt 1> /dev/null 2>&1; then
            cp "$sample_dir"/*.txt "$target_dir/"
        fi
        
        # Copy pypolca results if they exist
        if [[ -d "$sample_dir/pypolca" ]]; then
            find "$sample_dir/pypolca" -name "*.vcf" -exec cp {} "$target_dir/" \;
            find "$sample_dir/pypolca" -name "*report" -exec cp {} "$target_dir/" \;
        fi
    fi
done


# Save QUAST results for best assemblies improved with pypolca
echo "Saving QUAST results for best assemblies improved with pypolca..."
for sample_dir in "$workflow_path/$results_directory_name/best_assembly_qc"/*; do
    if [[ -d "$sample_dir" ]]; then
        # Extract sample name from path (best_assembly_qc/048ds -> 048ds)
        sample_name=$(basename "$sample_dir")
        
        # Check if quast_pypolca exists
        quast_dir="$sample_dir/quast_pypolca"
        if [[ ! -d "$quast_dir" ]]; then
            continue
        fi
        
        # Create descriptive directory name
        target_dir="$where_to_save/$name/${sample_name}/best_assembly_info_and_QC/quast_pypolca"
        
        echo "  Copying $sample_name QUAST results to $target_dir"
        mkdir -p "$target_dir"
        
        # Copy report files
        if ls "$quast_dir"/*report.tsv 1> /dev/null 2>&1; then
            cp "$quast_dir"/*report.tsv "$target_dir/" 2>/dev/null
        fi
    fi
done


# Save BUSCO results for best assemblies
echo "Saving best assembly BUSCO summary files..."
for sample_dir in "$workflow_path/$results_directory_name/best_assembly_qc"/*; do
    if [[ -d "$sample_dir" ]]; then
        # Extract sample name
        sample_name=$(basename "$sample_dir")
        
        # Process both busco types
        for busco_type in busco_general busco_specific; do
            busco_parent_dir="$sample_dir/${busco_type}_pypolca"
            
            if [[ ! -d "$busco_parent_dir" ]]; then
                continue
            fi
            
            # BUSCO creates a subdirectory named after the assembly
            # Find the actual BUSCO output directory (should be only one subdirectory)
            busco_dir=$(find "$busco_parent_dir" -mindepth 1 -maxdepth 1 -type d -name "BUSCO_*" | head -1)
            
            if [[ -z "$busco_dir" || ! -d "$busco_dir" ]]; then
                echo "  Warning: No BUSCO output directory found in $busco_parent_dir"
                continue
            fi
            
            # Create descriptive directory name
            target_dir="$where_to_save/$name/${sample_name}/best_assembly_info_and_QC/${busco_type}_pypolca"
            mkdir -p "$target_dir"
            
            echo "  Copying $sample_name ${busco_type}_pypolca summary files..."
            
            # Copy the main short_summary file (specific version)
            if find "$busco_dir" -maxdepth 1 -name "short_summary.specific.*.txt" -type f -print -quit | grep -q .; then
                find "$busco_dir" -maxdepth 1 -name "short_summary.specific.*.txt" -type f -exec cp {} "$target_dir/" \;
                echo "    Copied best assembly BUSCO summary"
            fi

            # Also copy full_table.tsv for downstream analysis
            if find "$busco_dir" -path "*/run_*/full_table.tsv" -type f -print -quit | grep -q .; then
                find "$busco_dir" -path "*/run_*/full_table.tsv" -type f -exec cp {} "$target_dir/${sample_name}_best_assembly_full_table.tsv" \;
                echo "    Copied ${busco_type} BUSCO full table"
            fi
        done
    fi
done


# Save MerquryFK results for best assemblies
echo "Saving best assembly MerquryFK results..."
for sample_dir in "$workflow_path/$results_directory_name/best_assembly_qc"/*; do
    if [[ -d "$sample_dir" ]]; then
        # Extract sample name  
        sample_name=$(basename "$sample_dir")
        
        merqury_dir="$sample_dir/merquryfk_pypolca"
        if [[ ! -d "$merqury_dir" ]]; then
            continue
        fi

        # Create target directory
        target_dir="$where_to_save/$name/${sample_name}/best_assembly_info_and_QC/merquryfk_pypolca"
        mkdir -p "$target_dir"
        
        echo "  Copying $sample_name best assembly MerquryFK results..."
        
        # Copy key MerquryFK files
        find "$merqury_dir" -maxdepth 1 -name "*.qv" -exec cp {} "$target_dir/" \;
        find "$merqury_dir" -maxdepth 1 -name "merquryfk.completeness.stats" -exec cp {} "$target_dir/" \;
        find "$merqury_dir" -maxdepth 1 -name "*ln.png" -exec cp {} "$target_dir/" \;
        find "$merqury_dir" -maxdepth 1 -name "*.bed" -exec cp {} "$target_dir/" \;
        
        echo "    Copied best assembly MerquryFK files"
    fi
done

# Save alignments of reads to best assemblies improved with pypolca
echo "Saving alignments and stats results for best assemblies improved with pypolca..."
for sample_dir in "$workflow_path/$results_directory_name/best_assembly_qc"/*; do
    if [[ -d "$sample_dir" ]]; then
        # Extract sample name from path
        sample_name=$(basename "$sample_dir")
        
        samtools_dir="$sample_dir/samtools_pypolca"
        if [[ ! -d "$samtools_dir" ]]; then
            continue
        fi
        
        # Create descriptive directory name
        target_dir="$where_to_save/$name/${sample_name}/best_assembly_info_and_QC/samtools_pypolca"

        echo "  Copying $sample_name samtools_pypolca results to $target_dir"
        mkdir -p "$target_dir"
        
        # Copy BAM files and indexes
        find "$samtools_dir" -name "*sorted.bam*" -exec cp {} "$target_dir/" \;
        # Copy stats files
        find "$samtools_dir" -name "*.txt" -exec cp {} "$target_dir/" \;
    fi
done

# Save coverage visualisation results for best assemblies improved with pypolca
echo "Saving coverage visualisation results for best assemblies improved with pypolca..."
for sample_dir in "$workflow_path/$results_directory_name/best_assembly_qc"/*; do
    if [[ -d "$sample_dir" ]]; then
        # Extract sample name from path
        sample_name=$(basename "$sample_dir")
        
        coverage_dir="$sample_dir/coverage_viz_pypolca"
        if [[ ! -d "$coverage_dir" ]]; then
            continue
        fi
        
        # Create descriptive directory name
        target_dir="$where_to_save/$name/${sample_name}/best_assembly_info_and_QC/coverage_viz_pypolca"

        echo "  Copying $sample_name coverage_viz_pypolca results to $target_dir"
        mkdir -p "$target_dir"
        
        # Copy PNG plots
        find "$coverage_dir" -name "*.png" -exec cp {} "$target_dir/" \;
        # Copy text files
        find "$coverage_dir" -name "*.txt" -exec cp {} "$target_dir/" \;
    fi
done



echo ""
echo "=============================================================="
echo "Copying final best assemblies to dedicated directory..."
echo "=============================================================="
echo ""

# Determine where to save final best assemblies
if [[ -z "$final_best_dir" ]]; then
    # If not specified, save in the run directory
    final_best_dir="$where_to_save/$name/FINAL_BEST_ASSEMBLIES"
else
    # Convert to absolute path if relative path provided
    final_best_dir=$(realpath -m "$final_best_dir")
fi

echo "Final best assemblies location: $final_best_dir"
mkdir -p "$final_best_dir"

# Copy best assembly and associated files for each sample
for sample_dir in "$workflow_path/$results_directory_name/best_assembly"/*; do
    if [[ -d "$sample_dir" ]]; then
        sample_name=$(basename "$sample_dir")
        
        echo "Processing final best assembly for sample: $sample_name"
        
        # Create sample directory in final location
        final_sample_dir="$final_best_dir/$sample_name"
        mkdir -p "$final_sample_dir"
        
        # Copy the best assembly fasta file
        best_assembly_file="$workflow_path/$results_directory_name/assemblies/$sample_name/${sample_name}_best_assembly_pypolca.fa"
        if [[ -f "$best_assembly_file" ]]; then
            cp "$best_assembly_file" "$final_sample_dir/${sample_name}_best_assembly.fa"
            echo "  ✓ Copied best assembly: ${sample_name}_best_assembly.fa"
        else
            echo "  🍄 Warning: Best assembly file not found: $best_assembly_file"
        fi
        
        # Copy the entire best_assembly_info_and_QC directory that was already saved
        qc_source="$where_to_save/$name/${sample_name}/best_assembly_info_and_QC"
        if [[ -d "$qc_source" ]]; then
            cp -r "$qc_source" "$final_sample_dir/"
            echo "  ✓ Copied QC results"
        else
            echo "  🍄 Warning: QC directory not found: $qc_source"
        fi
        
        # Copy logs for this sample and zip them
        logs_source="$where_to_save/$name/${name}_logs/$sample_name"
        if [[ -d "$logs_source" ]]; then
            cp -r "$logs_source" "$final_sample_dir/logs"
            # Create zip file in parent directory, then remove the logs folder
            (cd "$final_sample_dir" && zip -r -q logs.zip logs/)
            rm -rf "$final_sample_dir/logs"
            echo "  ✓ Copied and compressed logs"
        else
            echo "  🍄 Warning: Logs directory not found: $logs_source"
        fi

        # Copy the best assembly info file and rename for extract_metrics.py compatibility
        best_info_file="$final_sample_dir/best_assembly_info_and_QC/info/best_assembly.txt"
        if [[ -f "$best_info_file" ]]; then
            # copy as best_assembly_source.txt for backward compatibility with extract_metrics.py
            cp "$best_info_file" "$final_sample_dir/best_assembly_source.txt"
            winner=$(cat "$best_info_file")
            echo "  ✓ Best assembly was: $winner"
        fi
        
        echo ""
    fi
done

echo "=============================================================="
echo "Compressing logs and benchmarks..."
echo "=============================================================="
echo ""

# Zip the main logs directory
if [[ -d "$where_to_save/$name/${name}_logs" ]]; then
    echo "Compressing logs directory..."
    (cd "$where_to_save/$name" && zip -r -q "${name}_logs.zip" "${name}_logs/")
    rm -rf "$where_to_save/$name/${name}_logs"
    echo "  ✓ Logs compressed to ${name}_logs.zip"
fi

# Zip the benchmark directory
if [[ -d "$where_to_save/$name/${name}_benchmark" ]]; then
    echo "Compressing benchmark directory..."
    (cd "$where_to_save/$name" && zip -r -q "${name}_benchmark.zip" "${name}_benchmark/")
    rm -rf "$where_to_save/$name/${name}_benchmark"
    echo "  ✓ Benchmarks compressed to ${name}_benchmark.zip"
fi

echo ""
echo "=============================================================="
echo "Results saved successfully!"
echo "=============================================================="
echo ""
echo "Main results directory: $where_to_save/$name"
echo "Final best assemblies: $final_best_dir"
echo ""
echo "Summary of best assemblies:"
for sample_dir in "$final_best_dir"/*; do
    if [[ -d "$sample_dir" ]]; then
        sample_name=$(basename "$sample_dir")
        if [[ -f "$sample_dir/best_assembly_info_and_QC/info/best_assembly.txt" ]]; then
            winner=$(cat "$sample_dir/best_assembly_info_and_QC/info/best_assembly.txt")
            echo "  $sample_name: $winner"
        fi
    fi
done
echo ""

echo ""
echo "=============================================================="
echo "All done! Results are ready."
echo "Script completed at: $(date)"
echo "=============================================================="
echo ""
echo "Log file saved to: $log_file"
echo ""