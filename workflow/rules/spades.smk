# This rule runs spades assembler with the parameters specified in the config file.
import glob
import os


rule spades:
    input:
        forward_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R1.fq.gz",
        reverse_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R2.fq.gz",
        kmergenie_result=get_kmergenie_dependency,
    output:
        result_dir=directory(f"{output_dir}" + "{sample}/spades"),
        scaffolds=f"{output_dir}" + "{sample}/spades/scaffolds.fasta",
        link_assembly=f"{output_dir}" + "{sample}/assemblies/{sample}_megahit.fa",
    params:
        k=lambda wildcards: get_kmer_list(wildcards, "spades", "k"),
        optional_params=" ".join(
            f"{k} {v}" if v is not True else k
            for k, v in config["spades"]["optional_params"].items()
            if v and v is not False and v != ""
        ),
    threads: get_scaled_threads  # Use scaling function
    log:
        "logs/{sample}/spades.log",
    benchmark:
        "benchmark/{sample}/spades.txt"
    resources:
        mem_mb=get_scaled_mem,  # Use scaling function
    conda:
        "../envs/spades.yaml"
    shell:
        """
        echo "Running spades with the following command:" >> {log} 2>&1
        echo "spades.py -t {threads} -1 {input.forward_in} -2 {input.reverse_in} -o {output.result_dir} -k {params.k} {params.optional_params}" >> {log} 2>&1
        spades.py -t {threads} -1 {input.forward_in} -2 {input.reverse_in} -o {output.result_dir} -k {params.k} {params.optional_params} >> {log} 2>&1
        ln -srn {output.scaffolds} {output.link_assembly}
        """
