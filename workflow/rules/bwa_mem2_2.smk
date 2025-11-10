# This rule indexes the assemblies produced by the workflow, alignes the raw reads to the genome, and produces statistics
import glob
import os


rule bwa_mem2_2:
    input:
        assembly=f"{output_dir}" + "{sample}/assemblies/{sample}_best_assembly_pilon.fa",
        forward_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R1.fq.gz",
        reverse_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R2.fq.gz",
    output:
        sam=f"{output_dir}"
        + "{sample}/best_assembly_qc/bwa_mem2_pilon/{sample}_best_assembly_pilon.sam",
    threads: get_medium_high_threads
    resources:
        mem_mb=get_medium_high_mem,
        partition=config["medium_high"]["partition"],
    log:
        "logs/{sample}/bwa_mem2_best_assembly_pilon.log",
    benchmark:
        "benchmark/{sample}/bwa_mem2_best_assembly_pilon.txt"
    conda:
        "../envs/bwa_mem2.yaml"
    shell:
        """
        echo "Running bwa-mem2 index with the following command:" >> {log} 2>&1
        echo "bwa-mem2 index {input.assembly}" >> {log} 2>&1
        bwa-mem2 index {input.assembly} >> {log} 2>&1

        echo "Running bwa-mem2 mem with the following command:" >> {log} 2>&1
        echo "bwa-mem2 mem -t {threads} {input.assembly} {input.forward_in} {input.reverse_in} > {output.sam}" >> {log} 2>&1
        (bwa-mem2 mem -t {threads} {input.assembly} {input.forward_in} {input.reverse_in} > {output.sam}) >> {log} 2>&1
        """
