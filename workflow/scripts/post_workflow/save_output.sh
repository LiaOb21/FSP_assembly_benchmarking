#!/bin/bash

# This script is designed to save the output of the FSP assembly benchmarking pipeline easily.
# We will save: assembly fasta files, benchmarks, logs, QC results

# Function to show help
show_help() {
    cat << EOF
Usage: $0 -i <path_to_workflow_dir> -r <results_directory_name> -w <where_to_save> -b <batch_number> -t <reads_type> -k <kmer_strategy>

This script saves the output of the FSP assembly benchmarking pipeline.

Options:
  -i <path>   Path to the workflow root directory (FSP_assembly_benchmarking)
  -r <name>   Name of the directory containing the results of the workflow
  -w <path>   Path where the results should be saved
  -b <batch>  Batch number (e.g., EG1, batch2, B01, etc.). If you used multiple batches, use all the identifiers (e.g., EG1_EG2). If you don't have batches, use whatever you like in this field.
  -t <type>   Reads used: R1R2 or merged
  -k <kmer>   Kmer strategy: manual, kmergenie, or reads_length
  -h          Show this help message

The run name will be automatically generated as: <batch>_<reads_type>_<kmer_strategy>

Examples:
  $0 -i /path/to/FSP_assembly_benchmarking -r results -w /storage/final -b batch1 -t R1R2 -k manual
  → Creates: batch1_R1R2_manual

  $0 -i . -r results -w ./saved_output -b B02 -t merged -k kmergenie  
  → Creates: B02_merged_kmergenie

EOF
}

# Initialise variables
workflow_path=""
results_directory_name=""
where_to_save=""
batch_number=""
reads_type=""
kmer_strategy=""

# getopts to set user-defined variables
while getopts "i:r:w:b:t:k:h" opt; do
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
if [[ -z "$workflow_path" || -z "$results_directory_name" || -z "$where_to_save" || -z "$batch_number" || -z "$reads_type" || -z "$kmer_strategy" ]]; then
    echo "Error: All parameters are required!" >&2
    echo ""
    show_help
    exit 1
fi

# Generate the run name automatically
name="${batch_number}_${reads_type}_${kmer_strategy}"
echo ""
echo "Generated run name: $name"
echo ""

# Produce a simple text file containing the list of the samples used in the run
echo "Saving sample list..."
mkdir -p "$where_to_save/$name"
ls "$workflow_path/$results_directory_name" > "$where_to_save/$name/${name}_sample_list.txt"

# save fasta files. All the fasta files are located in $results_directory_name/$sample_name/assemblies (including the best assembly improved with pilon)
mkdir -p "$where_to_save/$name"

for sample_dir in "$workflow_path/$results_directory_name"/*/assemblies; do
    if [[ -d "$sample_dir" ]]; then
        # Extract sample name from path (e.g., results/048ds/assemblies -> 048ds)
        sample_name=$(basename "$(dirname "$sample_dir")")
        
        # Create target directory maintaining sample structure
        target_dir="$where_to_save/$name/$sample_name/assemblies"
        
        echo "  Copying $sample_name assemblies to $target_dir"
        mkdir -p "$target_dir"
        cp -rL "$sample_dir"/*fa "$target_dir/"
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
for busco_dir in "$workflow_path/$results_directory_name"/*/busco_*; do
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
                assembler_name=$(basename "$assembler_subdir")
                
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
                    cp $full_table_file "$target_dir/"
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
        cp -r "$sample_dir"/pilon/*changes "$target_dir/"
    fi
done


# Save QUAST results for each sample with sample name preserved
echo "Saving QUAST results for best assemblies improved with pilon..."
for sample_dir in "$workflow_path/$results_directory_name"/*/best_assembly_qc/quast_pilon; do
    if [[ -d "$sample_dir" ]]; then
        # Extract sample name from path
        sample_name=$(basename "$(dirname "$(dirname "$sample_dir")")")
        
        # Create descriptive directory name
        target_dir="$where_to_save/$name/${sample_name}/best_assembly_info_and_QC/quast_pilon"
        
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
            cp $full_table_file "$target_dir/"
            echo "    Copied $busco_type BUSCO full table"
        fi
    fi
done


# Save MerquryFK results - key files only
echo "Saving best assembly MerquryFK results..."
for merqury_dir in "$workflow_path/$results_directory_name"/*/best_assembly_qc/merquryfk_pilon; do
    if [[ -d "$merqury_dir" ]]; then
        # Extract sample name  
        sample_name=$(basename "$(dirname "$(dirname "$merqury_dir")")")

        # Create target directory
        target_dir="$where_to_save/$name/${sample_name}/best_assembly_info_and_QC/merquryfk_pilon"
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
        
        # Copy any summary files
        find "$merqury_dir" -maxdepth 1 -name "*.summary" -exec cp {} "$target_dir/" \;
        
        echo "    Copied best assembly MerquryFK files"
    fi
done

# Save alignments of reads to best assemblies improved with pilon (BAM files and indexes)
echo "Saving alignments and stats results for best assemblies improved with pilon..."
for sample_dir in "$workflow_path/$results_directory_name"/*/best_assembly_qc/bwa_mem2_samtools_pilon; do
    if [[ -d "$sample_dir" ]]; then
        # Extract sample name from path
        sample_name=$(basename "$(dirname "$(dirname "$sample_dir")")")
        
        # Create descriptive directory name
        target_dir="$where_to_save/$name/${sample_name}/best_assembly_info_and_QC/bwa_mem2_samtools_pilon"

        echo "  Copying $sample_name bwa_mem2_samtools_pilon results to $target_dir"
        mkdir -p "$target_dir"
        cp -r "$sample_dir"/* "$target_dir/"
    fi
done

# Save coverage visualisation results for best assemblies improved with pilon
echo "Saving coverage visualisation results for best assemblies improved with pilon..."
for sample_dir in "$workflow_path/$results_directory_name"/*/best_assembly_qc/coverage_viz_pilon; do
    if [[ -d "$sample_dir" ]]; then
        # Extract sample name from path
        sample_name=$(basename "$(dirname "$(dirname "$sample_dir")")")
        
        # Create descriptive directory name
        target_dir="$where_to_save/$name/${sample_name}/best_assembly_info_and_QC/coverage_viz_pilon"

        echo "  Copying $sample_name coverage_viz_pilon results to $target_dir"
        mkdir -p "$target_dir"
        cp -r "$sample_dir"/*png "$target_dir/"
    fi
done


echo "Results saved successfully to: $where_to_save/$name"