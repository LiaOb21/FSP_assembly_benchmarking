#!/bin/bash

# Select the best assembly per sample between multiple runs

# This script will iterate through multiple results directories, after saving them in the permanent directory.
# It will select the best assembly per sample based on QUAST and BUSCO scores, and copy the best assemblies to a final directory.

# Compare multiple directories with the same batch number and run identifier, but different reads types or kmer strategy.

# Function to show help
show_help() {
    cat << EOF
Usage: $0 -w <path_to_where_to_save> -r <final_best_assembly_results_directory_name> -b <batch_number> -n <run_id>

This script compares different runs best assemblies and saves the final output of the FSP assembly benchmarking pipeline to the final directory.
This script must be run after save_output.sh. The -w, -b, and -n options must match those used in save_output.sh to compare the correct runs where the same samples where assembled using different settings.

Options:
  -w <path>   Path where the results were previously saved using save_output.sh
  -r <name>   Name of the directory where the final best assemblies will be saved
  -b <batch>  Batch number (e.g., EG1, batch2, B01, etc.). This must be the same as the one used in save_output.sh
  -n <run_id> Run identifier (e.g., run1, subset_A, test, trial1, etc.). This must be the same as the one used in save_output.sh
  -h          Show this help message

Examples:
  $0 -w /path/to/FSP_assembly_benchmarking -r path/to/final_results -b batch1 -n run1
  → Compares: batch1_R1R2_manual_run1, batch1_merged_manual_run1, batch1_R1R2_kmergenie_run1, batch1_merged_kmergenie_run1, batch1_R1R2_reads_length_run1, batch1_merged_reads_length_run1
  → Saves the best assemblies per sample to: /path/to/FSP_assembly_benchmarking/path/to/final_results
EOF
}

# Parse command-line options
while getopts "w:r:b:n:h" opt; do
    case ${opt} in
        w ) where_to_save="$OPTARG" ;;
        r ) results_directory_name="$OPTARG" ;;
        b ) batch_number="$OPTARG" ;;
        n ) run_id="$OPTARG" ;;
        h ) show_help; exit 0 ;;
        \? ) echo "Invalid option: -$OPTARG" >&2; show_help; exit 1 ;;
        : ) echo "Option -$OPTARG requires an argument." >&2; show_help; exit 1 ;;
    esac
done

# Check if all required arguments are provided
if [[ -z "$where_to_save" || -z "$results_directory_name" || -z "$batch_number" || -z "$run_id" ]]; then
    echo "🍄 Error: Missing required arguments." >&2
    show_help
    exit 1
fi


# Function to copy assembly files and associated data
copy_assembly_files() {
    local best_run="$1"
    local best_assembler="$2"
    local sample="$3"
    local sample_output_dir="$4"
    
    # Copy the assembly
    local source_assembly="$where_to_save/$best_run/$sample/assemblies/${sample}_best_assembly_pypolca.fa"
    if [[ -f "$source_assembly" ]]; then
        cp "$source_assembly" "$sample_output_dir/${sample}_best_assembly.fa"
        echo "$best_run:$best_assembler" > "$sample_output_dir/best_assembly_source.txt"
        echo "    Assembly copied to: $sample_output_dir/${sample}_best_assembly.fa"
        
        # Copy the entire QC directory for the best assembly
        local source_qc_dir="$where_to_save/$best_run/$sample/best_assembly_info_and_QC"
        if [[ -d "$source_qc_dir" ]]; then
            cp -r "$source_qc_dir" "$sample_output_dir/best_assembly_info_and_QC"
            echo "    QC directory copied to: $sample_output_dir/best_assembly_info_and_QC"
        else
            echo "    🍄 Warning: QC directory not found: $source_qc_dir"
            echo "    Do you want to continue? (y/N)"
            read -r response
            if [[ ! "$response" =~ ^[Yy]$ ]]; then
                echo "Aborting."
                exit 1
            fi
        fi

        # Copy the logs for the best assembly
        local source_logs_dir="$where_to_save/$best_run/${best_run}_logs/$sample"
        if [[ -d "$source_logs_dir" ]]; then
            cp -r "$source_logs_dir" "$sample_output_dir/logs"
            echo "    Logs copied to: $sample_output_dir/logs"
        else
            echo "    🍄 Warning: Logs directory not found: $source_logs_dir"
            echo "    Do you want to continue? (y/N)"
            read -r response
            if [[ ! "$response" =~ ^[Yy]$ ]]; then
                echo "Aborting."
                exit 1
            fi
        fi

        # Copy the config file from the winning run
        local source_config="$where_to_save/$best_run/${best_run}_config.yml"
        if [[ -f "$source_config" ]]; then
            cp "$source_config" "$sample_output_dir/${sample}_best_assembly_config.yml"
            echo "    Config copied to: $sample_output_dir/${sample}_best_assembly_config.yml"
        else
            echo "    🍄 Warning: Config file not found: $source_config"
            echo "    Do you want to continue? (y/N)"
            read -r response
            if [[ ! "$response" =~ ^[Yy]$ ]]; then
                echo "Aborting."
                exit 1
            fi
        fi
    else
        echo "  🍄 Error: Source assembly not found: $source_assembly"
        echo "    Do you want to continue? (y/N)"
        read -r response
        if [[ ! "$response" =~ ^[Yy]$ ]]; then
            echo "Aborting."
            exit 1
        fi
    fi
}



# Each subdirectory in where_to_save matching is named as:
# <batch_number>_<reads_type>_<kmer_strategy>_<run_id>
# Find each subdirectory in the where_to_save directory that matches the batch number and run identifier
subdirs=()
while IFS= read -r -d $'\0' dir; do
    subdirs+=("$(basename "$dir")")
done < <(find "$where_to_save" -maxdepth 1 -type d -name "*${batch_number}*${run_id}*" -print0)

if [[ ${#subdirs[@]} -eq 0 ]]; then
    echo "🍄 Error: No matching subdirectories found for batch '$batch_number' and run ID '$run_id' in '$where_to_save'." >&2
    exit 1
else
    echo "✓ Found ${#subdirs[@]} matching directories:"
    for dir in "${subdirs[@]}"; do
        echo "  → $dir"
    done
    echo ""
fi

# Inside each subdirectory, check the *_sample_list.txt. These files must match between runs.
# Get the sample list from the first subdirectory
first_subdir="${subdirs[0]}"
first_sample_list_file="$where_to_save/$first_subdir/${first_subdir}_sample_list.txt"
if [[ ! -f "$first_sample_list_file" ]]; then
    echo "🍄 Error: Sample list file '$first_sample_list_file' not found."
    exit 1
fi
mapfile -t reference_samples < "$first_sample_list_file"

# Compare sample lists and check for missing samples
echo "Validating sample lists across runs..."
all_unique_samples=()
validation_warnings=0

# Collect all unique samples across all runs
for subdir in "${subdirs[@]}"; do
    sample_list_file="$where_to_save/$subdir/${subdir}_sample_list.txt"
    if [[ ! -f "$sample_list_file" ]]; then
        echo "🍄 Error: Sample list file '$sample_list_file' not found."
        exit 1
    fi
    
    mapfile -t temp_samples < "$sample_list_file"
    echo "  → $subdir: ${#temp_samples[@]} samples"
    
    # Add to unique samples list if not already present
    for sample in "${temp_samples[@]}"; do
        if [[ ! " ${all_unique_samples[*]} " =~ " $sample " ]]; then
            all_unique_samples+=("$sample")
        fi
    done
done

echo ""

# Check each sample against each run and report missing ones
for sample in "${all_unique_samples[@]}"; do
    for subdir in "${subdirs[@]}"; do
        sample_list_file="$where_to_save/$subdir/${subdir}_sample_list.txt"
        
        if ! grep -Fxq "$sample" "$sample_list_file"; then
            echo "🍄 Sample '$sample' is missing from run '$subdir'"
            ((validation_warnings++))
        fi
    done
done

# Ask user if they want to proceed when samples are missing
if [[ $validation_warnings -gt 0 ]]; then
    echo ""
    echo "🍄 Found $validation_warnings missing sample(s) across runs."
    echo "This might indicate failed analysis for some samples."
    echo "Do you want to proceed anyway? (y/N)"
    read -r response
    if [[ ! "$response" =~ ^[Yy]$ ]]; then
        echo "Operation cancelled. Please check your runs and try again."
        exit 1
    fi
    echo ""
    echo " Proceeding with missing samples - only available assemblies will be compared per sample."
else
    echo "✓ All samples are present in all runs!"
fi

echo ""
echo "✓ Sample validation complete! Will process ${#all_unique_samples[@]} samples:"
for i in "${!all_unique_samples[@]}"; do
    echo "  $((i+1)). ${all_unique_samples[$i]}"
done
echo ""

# Use all unique samples for processing
reference_samples=("${all_unique_samples[@]}")

# .............. to be continued ..............

# Set up final output directory
output_dir="$results_directory_name"
mkdir -p "$output_dir"

echo " Comparing best assemblies across runs..."
echo ""

# Process each sample
for sample in "${reference_samples[@]}"; do
    echo "Processing sample: $sample"
    
    sample_output_dir="$output_dir/$sample"
    # check if output directory for this sample exists, if not create it. If it exists, ask user before overwriting
    if [[ -d "$sample_output_dir" ]]; then
        echo "  🍄 Warning: Output directory for sample '$sample' already exists: $sample_output_dir"
        echo "  Do you want to overwrite it? (y/N)"
        read -r response
        if [[ ! "$response" =~ ^[Yy]$ ]]; then
            echo "  Skipping sample '$sample'."
            echo ""
            continue
        else
            echo "  Overwriting existing directory..."
            rm -rf "$sample_output_dir"
        fi
    fi

    if [[ ! -d "$sample_output_dir" ]]; then
        mkdir -p "$sample_output_dir"
    fi
    
    # Collect best assembly from each run that has this sample
    run_assemblies=()
    run_busco_scores=()
    
    for subdir in "${subdirs[@]}"; do
        sample_list_file="$where_to_save/$subdir/${subdir}_sample_list.txt"
        
        # Check if this run has this sample
        if grep -Fxq "$sample" "$sample_list_file"; then
            # Get the best assembly info for this sample in this run
            run_best_file="$where_to_save/$subdir/$sample/best_assembly_info_and_QC/info/best_assembly.txt"
            
            if [[ -f "$run_best_file" ]]; then
                best_assembler=$(cat "$run_best_file")
                echo "  → Run '$subdir': best assembler is '$best_assembler'"
                
                # Get BUSCO score for this best assembly
                busco_dir="$where_to_save/$subdir/$sample/best_assembly_info_and_QC/busco_specific_pypolca"
                busco_file=$(find "$busco_dir" -name "short_summary.specific.*.txt" 2>/dev/null | head -1)
                
                if [[ -f "$busco_file" ]]; then
                    busco_score=$(awk '/Complete and single-copy/ {print $1}' "$busco_file")
                    echo "    BUSCO score: $busco_score"
                    
                    # Store run info
                    run_assemblies+=("$subdir:$best_assembler")
                    run_busco_scores+=("$subdir:$best_assembler:$busco_score")
                else
                    echo "    🍄 Warning: No BUSCO file found for best assembly in run '$subdir'"
                    echo "    Do you want to continue? (y/N)"
                    read -r response
                    if [[ ! "$response" =~ ^[Yy]$ ]]; then
                        echo "Aborting."
                        exit 1
                    fi
                fi
            else
                echo "  🍄 Warning: No best_assembly.txt found in run '$subdir' for sample '$sample'"
                echo "  Do you want to continue? (y/N)"
                read -r response
                if [[ ! "$response" =~ ^[Yy]$ ]]; then
                    echo "Aborting."
                    exit 1
                fi
            fi
        fi
    done
    
    # Check if we have any assemblies to compare
    if [[ ${#run_assemblies[@]} -eq 0 ]]; then
        echo "  🍄 Error: No best assemblies found for sample '$sample' across any runs"
        continue
    fi
    
    if [[ ${#run_assemblies[@]} -eq 1 ]]; then
        # Only one run has this sample
        IFS=':' read -r best_run best_assembler <<< "${run_assemblies[0]}"
        echo "  ✓ Only one run has this sample - selecting '$best_assembler' from '$best_run'"
        
        copy_assembly_files "$best_run" "$best_assembler" "$sample" "$sample_output_dir"

    else
    
    # Multiple runs have this sample - compare BUSCO scores
    echo "  Comparing ${#run_assemblies[@]} assemblies across runs..."
    
    # Write BUSCO scores to temporary file for sorting
    temp_busco_file="$sample_output_dir/temp_busco_scores.txt"
    > "$temp_busco_file"
    
    for entry in "${run_busco_scores[@]}"; do
        IFS=':' read -r run assembler busco <<< "$entry"
        echo "$run $assembler $busco" >> "$temp_busco_file"
    done
    
    # Sort by BUSCO score (descending) and get the highest scoring assemblies
    sort -k3,3nr "$temp_busco_file" > "$sample_output_dir/sorted_busco_scores.txt"
    
    # Get the maximum BUSCO score
    max_busco=$(head -1 "$sample_output_dir/sorted_busco_scores.txt" | cut -d' ' -f3)
    
    # Get all assemblies with the maximum BUSCO score
    awk -v max="$max_busco" '$3 == max' "$sample_output_dir/sorted_busco_scores.txt" > "$sample_output_dir/best_busco_assemblies.txt"
    
    best_busco_count=$(wc -l < "$sample_output_dir/best_busco_assemblies.txt")
    
    if [[ $best_busco_count -eq 1 ]]; then
        # Single winner based on BUSCO
        read -r best_run best_assembler best_busco < "$sample_output_dir/best_busco_assemblies.txt"
        echo "  ✓ Best assembly: '$best_assembler' from '$best_run' (BUSCO: $best_busco)"
        
        copy_assembly_files "$best_run" "$best_assembler" "$sample" "$sample_output_dir"
    else
        # Tie in BUSCO scores - use auN as tiebreaker
        echo "  → Tie in BUSCO scores ($max_busco) - checking auN scores..."
        
        best_aun=0
        best_run=""
        best_assembler=""
        
        while read -r run assembler busco; do
            # Get auN score from QUAST report for this run
            quast_file="$where_to_save/$run/$sample/best_assembly_info_and_QC/quast_pypolca/report.tsv"
            if [[ -f "$quast_file" ]]; then
                # Since this is the best assembly QUAST report, there's only one column
                # Just extract the auN value directly
                aun_score=$(awk '/^auN/ {print $2}' "$quast_file")
                
                if [[ -n "$aun_score" ]] && awk "BEGIN {exit ($aun_score > $best_aun) ? 0 : 1}"; then
                    best_aun=$aun_score
                    best_run=$run
                    best_assembler=$assembler
                fi
                echo "    $run:$assembler - auN: $aun_score"
            else
                echo "    $run:$assembler - auN: not available (file not found)"
            fi
        done < "$sample_output_dir/best_busco_assemblies.txt"
        
        if [[ -n "$best_run" ]]; then
            echo "  ✓ Best assembly: '$best_assembler' from '$best_run' (BUSCO: $max_busco, auN: $best_aun)"
            
            copy_assembly_files "$best_run" "$best_assembler" "$sample" "$sample_output_dir"
        else
            echo "  🍄 Error: Could not determine best assembly for sample '$sample'"
            echo "    Do you want to continue? (y/N)"
            read -r response
            if [[ ! "$response" =~ ^[Yy]$ ]]; then
                echo "Aborting."
                exit 1
            fi
        fi
    fi
fi
    # Clean up temporary files
    rm -f "$sample_output_dir"/temp_*.txt
    echo ""
done

echo ""
echo "✓ Done! Final best assemblies saved in: $output_dir"
echo ""
echo "Summary of selected assemblies:"
for sample in "${reference_samples[@]}"; do
    if [[ -f "$output_dir/$sample/best_assembly_source.txt" ]]; then
        source_info=$(cat "$output_dir/$sample/best_assembly_source.txt")
        echo "  $sample: $source_info"
    fi
done