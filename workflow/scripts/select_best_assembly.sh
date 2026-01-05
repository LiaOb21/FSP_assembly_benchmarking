#!/bin/bash

# This script selects the best assembly based on BUSCO scores with auN as tiebreaker

# Function to show help
show_help() {
    cat << EOF
Usage: $0 -s <sample> -r <snakemake_results_directory> -o <output_directory>

This script finds the best assembly for each sample based on:
1. Highest complete and single-copy BUSCOs
2. If multiple assemblies have the highest BUSCOs, the best assembly is chosen among those with the best BUSCOs based on the auN value from QUAST

Options:
  -s <sample>  Sample name
  -r <path>   Path to snakemake results directory containing sample folders
  -o <path>   Output directory (optional, default: current directory, recommended: snakemake_results_directory/best_assemblies)
  -h          Show this help message

Example:
  $0 -r results
  $0 -s 048ds -r results -o results/best_assemblies_output 

EOF
}

# Initialize variables
sample=""
results_dir=""
output_dir="."

# Parse command line arguments
while getopts "s:r:o:h" opt; do
    case $opt in
        s) sample="$OPTARG"
           echo "Sample: $OPTARG"
        ;;
        r) results_dir="$OPTARG"
           echo "Results directory: $OPTARG"
        ;;
        o) output_dir="$OPTARG"
           echo "Output directory: $OPTARG"
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

# Check if results directory is provided
if [[ -z "$results_dir" ]]; then
    echo "Error: Results directory (-r) is required!" >&2
    show_help
    exit 1
fi

# Check if results directory exists
if [[ ! -d "$results_dir" ]]; then
    echo "Error: Results directory '$results_dir' does not exist!" >&2
    exit 1
fi

# Handle existing output directory - overwrite if it exists
if [[ -d "$output_dir" ]]; then
    echo "Output directory '$output_dir' already exists. Overwriting..."
    rm -rf "$output_dir"
fi

# Create output directory
mkdir -p "$output_dir"

echo "processing sample: $sample"

echo "=== Searching for assemblies across all strategies ==="

# Extract BUSCO complete and single-copy for each assembly across all strategies and read types
for reads_type_dir in "${results_dir}"*/; do
    if [[ -d "$reads_type_dir" ]]; then
        reads_type=$(basename "$reads_type_dir")
        
        # Skip non-reads-type directories (like 'assemblies', 'best_assembly', etc.)
        if [[ "$reads_type" != "R1R2" && "$reads_type" != "merged" ]]; then
            continue
        fi
        
        echo "  Checking reads_type: $reads_type" >&2
        
        for strategy_dir in "$reads_type_dir"/*/; do
            if [[ -d "$strategy_dir" ]]; then
                strategy=$(basename "$strategy_dir")
                sample_results_dir="$strategy_dir/$sample"
                
                echo "    Checking strategy: $strategy" >&2
                
                for assembler_dir in "$sample_results_dir"/busco_specific/*/; do
                    if [[ -d "$assembler_dir" ]]; then
                        assembler=$(basename "$assembler_dir")
                        for file in "$assembler_dir"/BUSCO_*/short_summary.specific.*.txt; do
                            if [[ -f "$file" ]]; then
                                # Include reads_type and strategy in output
                                busco_score=$(awk '/Complete and single-copy/ {print $1}' "$file")
                                echo "      Found: ${reads_type}_${strategy}_${assembler} (BUSCO: $busco_score)" >&2
                                awk -v assembler="$assembler" -v strategy="$strategy" -v reads_type="$reads_type" \
                                    '/Complete and single-copy/ {print reads_type"_"strategy"_"assembler, $1}' "$file" \
                                    >> "$output_dir/complete_single_copy_buscos.txt"
                            fi
                        done
                    fi
                done
            fi
        done
    fi
done

if [[ ! -f "$output_dir/complete_single_copy_buscos.txt" ]]; then
    echo "  No BUSCO data found for sample $sample"
    exit 1
fi

echo ""
echo "=== All assemblies with BUSCO scores ==="
cat "$output_dir/complete_single_copy_buscos.txt"
echo ""

# Sort the results by the number of complete and single-copy BUSCOs in descending order
# extract max value from the first line
# print all lines with that max value to best_buscos.txt
sort -k2,2nr "$output_dir/complete_single_copy_buscos.txt" | awk 'NR==1{max=$2} $2==max' > "$output_dir/best_buscos.txt"

echo "=== Assemblies with highest BUSCO scores ==="
cat "$output_dir/best_buscos.txt"
echo ""

# Check if only one assembly has the highest BUSCO score
best_busco_count=$(wc -l < "$output_dir/best_buscos.txt")

if [[ $best_busco_count -eq 1 ]]; then
    # Single winner - print the assembly name
    best_assembly=$(cut -d' ' -f1 "$output_dir/best_buscos.txt")
    echo "$best_assembly is the best assembly for sample $sample based on BUSCO score."
    echo "linking ${results_dir}assemblies/$sample/${sample}_${best_assembly}.fa to $output_dir"
    ln -srn "${results_dir}assemblies/$sample/${sample}_${best_assembly}.fa" "$output_dir/${sample}_best_assembly.fa"
    echo "$best_assembly" > "$output_dir/best_assembly.txt"
    exit 0
fi

# Extract Assembly and aUN values from quast report for this sample and transpose
echo "  Multiple assemblies with the highest BUSCOs, checking auN..."
echo "Extracting auN scores..."

# Extract auN from QUAST report
quast_file="${results_dir}/quast/$sample/report.txt"

if [[ -f "$quast_file" ]]; then
    echo "  Found QUAST report at: $quast_file" >&2
    # Extract Assembly and auN values from quast report
    # Assembly names format: sample_reads_type_strategy_assembler
    awk '
        /Assembly/ {
            for(i=2;i<=NF;i++) {
                # Remove .fa extension
                gsub(/\.fa$/, "", $i);
                
                # Extract everything starting from R1R2_ or merged_
                # This handles any sample name format including underscores
                if (match($i, /(R1R2|merged)_[^_]+_[^_]+$/)) {
                    extracted = substr($i, RSTART);
                    a[i-1] = extracted;
                }
            }
        }
        /auN/ {
            for(i=2;i<=NF;i++) {
                if (a[i-1] != "") {
                    print a[i-1], $i
                }
            }
        }
    ' "$quast_file" >> "$output_dir/auN_quast.txt"
fi

# Get auN scores only for the tied BUSCO assemblers
# Get the list of assemblers with the highest BUSCOs
tied_assemblers=$(cut -d' ' -f1 "$output_dir/best_buscos.txt")
# set up variables to track best auN and corresponding assembly
best_aun=0
best_assembly=""

# for each assembler in tied_assemblers list, get its auN score and compare
# if higher than current best_aun, update best_aun and best_assembly
while read assembler; do
    aun_score=$(grep "^$assembler " "$output_dir/auN_quast.txt" | cut -d' ' -f2)
    if [[ -n "$aun_score" ]] && awk "BEGIN {exit ($aun_score > $best_aun) ? 0 : 1}"; then
        best_aun=$aun_score
        best_assembly=$assembler
    fi
done <<< "$tied_assemblers"

echo ""
echo "=== auN scores for tied assemblies ==="
grep -F -f <(cut -d' ' -f1 "$output_dir/best_buscos.txt") "$output_dir/auN_quast.txt"
echo ""
    
if [[ -n "$best_assembly" ]]; then
    echo "$best_assembly is the best assembly for sample $sample based on auN score." >&2
    echo "  Best auN: $best_aun" >&2
    
    echo "linking ${results_dir}assemblies/$sample/${sample}_${best_assembly}.fa to $output_dir" >&2
    ln -srn "${results_dir}assemblies/$sample/${sample}_${best_assembly}.fa" "$output_dir/${sample}_best_assembly.fa"
    echo "$best_assembly" > "$output_dir/best_assembly.txt"
    exit 0
fi

echo "Done! Results saved in $output_dir"
