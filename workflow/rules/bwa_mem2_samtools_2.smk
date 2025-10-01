# This rule indexes the assemblies produced by the workflow, alignes the raw reads to the genome, and produces statistics
import glob
import os


rule bwa_mem2_samtools_2:
    input:
        assembly=f"{output_dir}" + "{sample}/assemblies/{sample}_best_assembly_pilon.fa",
        forward_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R1.fq.gz",
        reverse_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R2.fq.gz",
    output:
        sam=f"{output_dir}"
        + "{sample}/best_assembly_qc/bwa_mem2_samtools_pilon/{sample}_best_assembly_pilon.sam",
        bam=f"{output_dir}"
        + "{sample}/best_assembly_qc/bwa_mem2_samtools_pilon/{sample}_best_assembly_pilon.bam",
        sorted_bam=f"{output_dir}"
        + "{sample}/best_assembly_qc/bwa_mem2_samtools_pilon/{sample}_best_assembly_pilon_sorted.bam",
        coverage_stats=f"{output_dir}"
        + "{sample}/best_assembly_qc/bwa_mem2_samtools_pilon/{sample}_best_assembly_pilon_coverage_stats.txt",
        flagstat=f"{output_dir}"
        + "{sample}/best_assembly_qc/bwa_mem2_samtools_pilon/{sample}_best_assembly_pilon_flagstat.txt",
    threads: get_high_threads
    resources:
        mem_mb=get_high_mem,
        partition=config["high"]["partition"],
    log:
        "logs/{sample}/bwa_mem2_best_assembly_pilon.log",
    benchmark:
        "benchmark/{sample}/bwa_mem2_best_assembly_pilon.txt"
    conda:
        "../envs/bwa_mem2_samtools.yaml"
    shell:
        """
        echo "Running bwa-mem2 index with the following command:" >> {log} 2>&1
        echo "bwa-mem2 index {input.assembly}" >> {log} 2>&1
        bwa-mem2 index {input.assembly} >> {log} 2>&1

        echo "Running bwa-mem2 mem with the following command:" >> {log} 2>&1
        echo "bwa-mem2 mem -t {threads} {input.assembly} {input.forward_in} {input.reverse_in} > {output.sam}" >> {log} 2>&1
        (bwa-mem2 mem -t {threads} {input.assembly} {input.forward_in} {input.reverse_in} > {output.sam}) >> {log} 2>&1

        echo "Running samtools view with the following command:" >> {log} 2>&1
        echo "samtools view -bS {output.sam} -@ {threads}  > {output.bam}" >> {log} 2>&1
        (samtools view -bS {output.sam} -@ {threads} > {output.bam}) >> {log} 2>&1

        echo "Running samtools sort with the following command:" >> {log} 2>&1
        echo "samtools sort {output.bam} -@ {threads} > {output.sorted_bam}" >> {log} 2>&1
        (samtools sort {output.bam} -@ {threads} > {output.sorted_bam}) >> {log} 2>&1

        echo "Indexing the sorted bam file with samtools index" >> {log} 2>&1
        samtools index {output.sorted_bam} 

        echo "Calculating coverage statistics with samtools coverage:" >> {log} 2>&1
        (samtools coverage {output.sorted_bam} > {output.coverage_stats}) 2>> {log}

        echo "Calculating mapping statistics with samtools flagstats:" >> {log} 2>&1
        (samtools flagstat {output.sorted_bam} > {output.flagstat}) 2>> {log}        
        """
