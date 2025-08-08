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
        expand(f"{output_dir}" + "{sample}/spades/scaffolds.fasta", sample=SAMPLES),
        expand(f"{output_dir}" + "{sample}/megahit/", sample=SAMPLES),
        expand(f"{output_dir}" + "{sample}/abyss/abyss-scaffolds.fa", sample=SAMPLES),
        expand(f"{output_dir}" + "{sample}/sparseassembler/SuperContigs.txt", sample=SAMPLES),
        expand(f"{output_dir}" + "{sample}/minia/{sample}.contigs.fa", sample=SAMPLES),
        expand(f"{output_dir}" + "{sample}/masurca/masurca_config.txt", sample=SAMPLES),
        expand(f"{output_dir}" + "{sample}/masurca/CA/primary.genome.scf.fasta", sample=SAMPLES),
        expand(f"{output_dir}" + "{sample}/busco_general/{assembler}", sample=SAMPLES, assembler=ASSEMBLERS),
        expand(f"{output_dir}" + "{sample}/busco_specific/{assembler}", sample=SAMPLES, assembler=ASSEMBLERS),
        expand(f"{output_dir}" + "{sample}/quast", sample=SAMPLES),
        expand(f"{output_dir}" + "{sample}/fastk/fastk_table.ktab", sample=SAMPLES),
        expand(f"{output_dir}" + "{sample}/merquryfk/{assembler}/merquryfk.completeness.stats", sample=SAMPLES, assembler=ASSEMBLERS),
        expand(f"{output_dir}" + "{sample}/merquryfk/{assembler}/merquryfk.qv", sample=SAMPLES, assembler=ASSEMBLERS),
        expand(f"{output_dir}" + "{sample}/bwa_mem2_samtools/{assembler}/{sample}_{assembler}_sorted.bam", sample=SAMPLES, assembler=ASSEMBLERS),
        expand(f"{output_dir}" + "{sample}/coverage_viz/{assembler}/{sample}_{assembler}_coverage_plot.png", sample=SAMPLES, assembler=ASSEMBLERS),
        expand(f"{output_dir}" + "{sample}/pilon/{assembler}/{sample}_{assembler}_pilon.fasta", sample=SAMPLES, assembler=ASSEMBLERS),
        expand(f"{output_dir}" + "{sample}/busco_general_pilon/{assembler}", sample=SAMPLES, assembler=ASSEMBLERS),
        expand(f"{output_dir}" + "{sample}/busco_specific_pilon/{assembler}", sample=SAMPLES, assembler=ASSEMBLERS),
        expand(f"{output_dir}" + "{sample}/merquryfk_pilon/{assembler}/merquryfk.completeness.stats", sample=SAMPLES, assembler=ASSEMBLERS),
        expand(f"{output_dir}" + "{sample}/merquryfk_pilon/{assembler}/merquryfk.qv", sample=SAMPLES, assembler=ASSEMBLERS),
        expand(f"{output_dir}" + "{sample}/bwa_mem2_samtools_pilon/{assembler}/{sample}_{assembler}_sorted.bam", sample=SAMPLES, assembler=ASSEMBLERS),
        expand(f"{output_dir}" + "{sample}/coverage_viz_pilon/{assembler}/{sample}_{assembler}_coverage_plot.png", sample=SAMPLES, assembler=ASSEMBLERS),
        expand(f"{output_dir}" + "{sample}/mapDamage2/{assembler}", sample=SAMPLES, assembler=ASSEMBLERS),
        expand(f"{output_dir}" + "{sample}/mapDamage2_pilon/{assembler}", sample=SAMPLES, assembler=ASSEMBLERS),
    ]
    
    # Add kmergenie outputs only if mode is "auto"
    if config["kmer_strategy"]["mode"] == "auto":
        inputs.extend([
            expand(f"{output_dir}" + "{sample}/kmergenie/{sample}_report.html", sample=SAMPLES),
            expand(f"{output_dir}" + "{sample}/kmergenie/{sample}_best_kmer.txt", sample=SAMPLES),
        ])
    
    return inputs


# Function to read kmer size for the assemblers dynamic based on config.yml
# If auto mode, read from kmergenie output.
# kmer predicted by kmergenie must be a value between 15 and 127
# If manual mode, use config values.

def get_kmer_list(wildcards, assembler, config_key):
    """
    Get k-mer list for assemblers that use multiple k-mer sizes.
    For spades and megahit.
    """
    # Default k-mer list (same for both spades and megahit)
    default_k_list = [21, 27, 29, 31, 33, 35, 37, 41, 51, 61]
    
    if config["kmer_strategy"]["mode"] == "auto":
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
    if config["kmer_strategy"]["mode"] == "auto":
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


# This function is used to add kmergenie as dependency when kmer_strategy is auto (i.e. we want to run kmergenie)
def get_kmergenie_dependency(wildcards):
    """Return kmergenie dependency if in auto mode, empty list otherwise"""
    if config["kmer_strategy"]["mode"] == "auto":
        return f"{output_dir}{wildcards.sample}/kmergenie/{wildcards.sample}_best_kmer.txt"
    else:
        return []