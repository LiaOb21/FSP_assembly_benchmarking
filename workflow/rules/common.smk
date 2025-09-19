# resource scaling functions

def get_scaled_mem(wildcards, attempt):
    base_mem = config["mem_mb"]
    return min(base_mem * attempt, 250000)  # Cap at 250GB


def get_scaled_threads(wildcards, attempt):
    base_threads = config["threads"]
    return min(
        base_threads + (attempt - 1) * 4, 64
    )  # Add 4 threads per retry, cap at 64

# Make inputs for rule all dynamic based on k-mer strategy

def get_all_inputs():
    inputs = [
#        expand(f"{output_dir}" + "{sample}/spades/scaffolds.fasta", sample=SAMPLES),
#        expand(f"{output_dir}" + "{sample}/megahit/", sample=SAMPLES),
#        expand(f"{output_dir}" + "{sample}/abyss/abyss-scaffolds.fa", sample=SAMPLES),
#        expand(f"{output_dir}" + "{sample}/sparseassembler/SuperContigs.txt", sample=SAMPLES),
#        expand(f"{output_dir}" + "{sample}/minia/{sample}.contigs.fa", sample=SAMPLES),
#        expand(f"{output_dir}" + "{sample}/masurca/masurca_config.txt", sample=SAMPLES),
#        expand(f"{output_dir}" + "{sample}/masurca/CA/primary.genome.scf.fasta", sample=SAMPLES),
#        expand(f"{output_dir}" + "{sample}/busco_general/{assembler}", sample=SAMPLES, assembler=ASSEMBLERS),
#        expand(f"{output_dir}" + "{sample}/busco_specific/{assembler}", sample=SAMPLES, assembler=ASSEMBLERS),
#        expand(f"{output_dir}" + "{sample}/quast", sample=SAMPLES),
#        expand(f"{output_dir}" + "{sample}/fastk/fastk_table.ktab", sample=SAMPLES),
        expand(f"{output_dir}" + "{sample}/merquryfk/{assembler}/merquryfk.completeness.stats", sample=SAMPLES, assembler=ASSEMBLERS),
        expand(f"{output_dir}" + "{sample}/merquryfk/{assembler}/merquryfk.qv", sample=SAMPLES, assembler=ASSEMBLERS),
#        expand(f"{output_dir}" + "{sample}/best_assembly/{sample}_best_assembly.fa", sample=SAMPLES),
#        expand(f"{output_dir}" + "{sample}/best_assembly/bwa_mem2_samtools/{sample}_best_assembly_sorted.bam", sample=SAMPLES),
#        expand(f"{output_dir}" + "{sample}/coverage_viz/{assembler}/{sample}_{assembler}_coverage_plot.png", sample=SAMPLES, assembler=ASSEMBLERS),
#        expand(f"{output_dir}" + "{sample}/best_assembly/pilon/{sample}_best_assembly_pilon.fasta", sample=SAMPLES),
        expand(f"{output_dir}" + "{sample}/best_assembly_qc/busco_general_pilon", sample=SAMPLES, assembler=ASSEMBLERS),
        expand(f"{output_dir}" + "{sample}/best_assembly_qc/busco_specific_pilon", sample=SAMPLES, assembler=ASSEMBLERS),
        expand(f"{output_dir}" + "{sample}/best_assembly_qc/merquryfk_pilon/merquryfk.completeness.stats", sample=SAMPLES, assembler=ASSEMBLERS),
        expand(f"{output_dir}" + "{sample}/best_assembly_qc/merquryfk_pilon/merquryfk.qv", sample=SAMPLES, assembler=ASSEMBLERS),
        expand(f"{output_dir}" + "{sample}/best_assembly_qc/quast_pilon", sample=SAMPLES),
        expand(f"{output_dir}" + "{sample}/best_assembly_qc/bwa_mem2_samtools_pilon/{sample}_best_assembly_pilon_sorted.bam", sample=SAMPLES, assembler=ASSEMBLERS),
        expand(f"{output_dir}" + "{sample}/best_assembly_qc/coverage_viz_pilon/{sample}_best_assembly_pilon_coverage_summary.txt", sample=SAMPLES),
    ]
    
    # Add kmergenie outputs only if mode is "kmergenie"
    if config["kmer_strategy"]["mode"] == "kmergenie":
        inputs.extend([
            expand(f"{output_dir}" + "{sample}/kmergenie/{sample}_report.html", sample=SAMPLES),
            expand(f"{output_dir}" + "{sample}/kmergenie/{sample}_best_kmer.txt", sample=SAMPLES),
        ])
    
    return inputs


# Function to read kmer size for the assemblers dynamic based on config.yml
# If kmergenie mode, read from kmergenie output.
# kmer predicted by kmergenie must be a value between 15 and 127
# If manual mode, use config values.

def get_kmer_list(wildcards, assembler, config_key):
    """
    Get k-mer list for assemblers that use multiple k-mer sizes.
    For spades and megahit.
    """
    # Default k-mer list (same for both spades and megahit)
    default_k_list = [21, 29, 33, 39, 55, 59, 79, 99, 119, 127]
    
    if config["kmer_strategy"]["mode"] == "kmergenie":
        # Read best k-mer from kmergenie output
        best_kmer_file = f"{output_dir}{wildcards.sample}/kmergenie/{wildcards.sample}_best_kmer.txt"
        with open(best_kmer_file, 'r') as f:
            best_k = int(f.read().strip())
        
        # Add best k if not present and within range (15-127)
        if best_k not in default_k_list and 15 <= best_k <= 127 and best_k % 2 == 1:
            k_list = default_k_list + [best_k]
            k_list.sort()
        else:
            # If kmergenie prediction is outside range, use default list
            k_list = default_k_list
        
        return ",".join(map(str, k_list))
    else:
        # Manual mode - use config value (user must set this)
        return config[assembler][config_key]


def get_single_kmer(wildcards, assembler, config_key):
    """
    Get single k-mer value for assemblers that use only one k-mer size.
    For abyss, sparseassembler, minia, masurca.
    """
    if config["kmer_strategy"]["mode"] == "kmergenie":
        # Read best k-mer from kmergenie output
        best_kmer_file = f"{output_dir}{wildcards.sample}/kmergenie/{wildcards.sample}_best_kmer.txt"
        with open(best_kmer_file, 'r') as f:
            best_k = int(f.read().strip())
        
        # Use kmergenie prediction if within range (15-127), otherwise use 25
        if 15 <= best_k <= 127 and best_k % 2 == 1:
            return str(best_k)
        else:
            return "25"  # Fallback to 25 if outside range
    else:
        # Manual mode - use config value (user must set this)
        return str(config[assembler][config_key])


# This function is used to add kmergenie as dependency when kmer_strategy is kmergenie (i.e. we want to run kmergenie)
def get_kmergenie_dependency(wildcards):
    """Return kmergenie dependency if in kmergenie mode, empty list otherwise"""
    if config["kmer_strategy"]["mode"] == "kmergenie":
        return f"{output_dir}{wildcards.sample}/kmergenie/{wildcards.sample}_best_kmer.txt"
    else:
        return []