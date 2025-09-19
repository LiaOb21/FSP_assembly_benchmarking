# This rule indexes the assemblies produced by the workflow, alignes the raw reads to the genome, and produces statistics
import glob
import os


rule bwa_mem2_samtools:
    input:
        assembly=f"{output_dir}" + "{sample}/best_assembly/{sample}_best_assembly.fa",
        forward_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R1.fq.gz",
        reverse_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R2.fq.gz",
    output:
        sam=f"{output_dir}" + "{sample}/best_assembly/bwa_mem2_samtools/{sample}_best_assembly.sam",
        bam=f"{output_dir}" + "{sample}/best_assembly/bwa_mem2_samtools/{sample}_best_assembly.bam",
        sorted_bam=f"{output_dir}" + "{sample}/best_assembly/bwa_mem2_samtools/{sample}_best_assembly_sorted.bam",
        bam_index=f"{output_dir}" + "{sample}/best_assembly/bwa_mem2_samtools/{sample}_best_assembly_sorted.bam.bai", 
    threads: get_scaled_threads  # Use scaling function
    log:
        "logs/{sample}/bwa_mem2_best_assembly.log",
    benchmark:
        "benchmark/{sample}/bwa_mem2_best_assembly.txt"
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

        echo "Indexing the sorted bam file with samtools index" >> {log} 2>&1
        samtools index {output.sorted_bam}      
        """
