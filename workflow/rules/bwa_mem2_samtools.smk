# This rule indexes the assemblies produced by the workflow
import glob
import os


rule bwa_mem2_samtools:
    input:
        assembly=f"{output_dir}" + "assemblies/{sample}/{sample}_{assembler}.fa",
        forward_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R1.fq.gz",
        reverse_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R2.fq.gz",
    output:
        sam=f"{output_dir}" + "{sample}/bwa_mem2_samtools/{assembler}/{sample}_{assembler}.sam",
        bam=f"{output_dir}" + "{sample}/bwa_mem2_samtools/{assembler}/{sample}_{assembler}.bam",
        sorted_bam=f"{output_dir}" + "{sample}/bwa_mem2_samtools/{assembler}/{sample}_{assembler}_sorted.bam",
        coverage_stats=f"{output_dir}" + "{sample}/bwa_mem2_samtools/{assembler}/{sample}_{assembler}_coverage_stats.txt",
        flagstat=f"{output_dir}" + "{sample}/bwa_mem2_samtools/{assembler}/{sample}_{assembler}_flagstat.txt",
    threads: get_scaled_threads  # Use scaling function
    log:
        "logs/{sample}/bwa_mem2_{assembler}.log",
    benchmark:
        "benchmark/{sample}/bwa_mem2_{assembler}.txt"
    resources:
        mem_mb=get_scaled_mem,  # Use scaling function
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

        echo "Calculating coverage statistics with samtools coverage:" >> {log} 2>&1
        (samtools coverage {output.sorted_bam} > {output.coverage_stats}) 2>> {log}

        echo "Calculating mapping statistics with samtools flagstats:" >> {log} 2>&1
        (samtools flagstat {output.sorted_bam} > {output.flagstat}) 2>> {log}        
        """
