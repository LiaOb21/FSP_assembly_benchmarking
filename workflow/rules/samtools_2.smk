# This rule indexes the assemblies produced by the workflow, alignes the raw reads to the genome, and produces statistics
import glob
import os


rule samtools_2:
    input:
        sam=f"{output_dir}"
        + "{sample}/best_assembly_qc/bwa_mem2_pilon/{sample}_best_assembly_pilon.sam",
    output:
        bam=f"{output_dir}"
        + "{sample}/best_assembly_qc/samtools_pilon/{sample}_best_assembly_pilon.bam",
        sorted_bam=f"{output_dir}"
        + "{sample}/best_assembly_qc/samtools_pilon/{sample}_best_assembly_pilon_sorted.bam",
        coverage_stats=f"{output_dir}"
        + "{sample}/best_assembly_qc/samtools_pilon/{sample}_best_assembly_pilon_coverage_stats.txt",
        flagstat=f"{output_dir}"
        + "{sample}/best_assembly_qc/samtools_pilon/{sample}_best_assembly_pilon_flagstat.txt",
    threads: get_medium_high_threads
    resources:
        mem_mb=get_medium_high_mem,
        partition=config["medium_high"]["partition"],
    log:
        "logs/{sample}/samtools_best_assembly_pilon.log",
    benchmark:
        "benchmark/{sample}/samtools_best_assembly_pilon.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        echo "Running samtools view with the following command:" >> {log} 2>&1
        echo "samtools view -b {input.sam} -@ {threads}  > {output.bam}" >> {log} 2>&1
        (samtools view -b {input.sam} -@ {threads} > {output.bam}) >> {log} 2>&1

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
