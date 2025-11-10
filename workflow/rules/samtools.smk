# This rule indexes the assemblies produced by the workflow, alignes the raw reads to the genome, and produces statistics
import glob
import os


rule samtools:
    input:
        sam=f"{output_dir}"
        + "{sample}/best_assembly/bwa_mem2/{sample}_best_assembly.sam",
    output:
        bam=f"{output_dir}"
        + "{sample}/best_assembly/samtools/{sample}_best_assembly.bam",
        name_sorted_bam=f"{output_dir}"
        + "{sample}/best_assembly/samtools/{sample}_best_assembly_name_sorted.bam",
        bam_fixmate=f"{output_dir}"
        + "{sample}/best_assembly/samtools/{sample}_best_assembly_fixmate.bam",
        sorted_bam=f"{output_dir}"
        + "{sample}/best_assembly/samtools/{sample}_best_assembly_fixmate_sorted.bam",
        bam_markdup=f"{output_dir}"
        + "{sample}/best_assembly/samtools/{sample}_best_assembly_fixmate_sorted_markdup.bam",
        bam_index=f"{output_dir}"
        + "{sample}/best_assembly/samtools/{sample}_best_assembly_fixmate_sorted_markdup.bam.bai",
    threads: get_medium_high_threads
    resources:
        mem_mb=get_medium_high_mem,
        partition=config["medium_high"]["partition"],
    log:
        "logs/{sample}/samtools_best_assembly.log",
    benchmark:
        "benchmark/{sample}/samtools_best_assembly.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        echo "Running samtools view with the following command:" >> {log} 2>&1
        echo "samtools view -b -F 4 -q 20 {input.sam} -@ {threads}  > {output.bam}" >> {log} 2>&1
        echo "Filtering unmapped reads (-F 4) and low quality alignments (-q 20) from SAM to BAM format" >> {log} 2>&1
        (samtools view -b -F 4 -q 20 {input.sam} -@ {threads} > {output.bam}) >> {log} 2>&1

        echo "Sorting bam file according to reads name:" >> {log} 2>&1
        echo "samtools sort -n {output.bam} -@ {threads} > {output.name_sorted_bam}" >> {log} 2>&1
        (samtools sort -n {output.bam} -@ {threads} > {output.name_sorted_bam}) >> {log} 2>&1

        echo "Running samtools fixmate for subsequent marking duplicates:" >> {log} 2>&1
        echo "samtools fixmate -m {output.name_sorted_bam} {output.bam_fixmate}" >> {log} 2>&1
        (samtools fixmate -m {output.name_sorted_bam} {output.bam_fixmate}) >> {log} 2>&1

        echo "Sorting fixmate bam file according to coordinates:" >> {log} 2>&1
        echo "samtools sort {output.bam_fixmate} -@ {threads} > {output.sorted_bam}" >> {log} 2>&1
        (samtools sort {output.bam_fixmate} -@ {threads} > {output.sorted_bam}) >> {log} 2>&1

        echo "Marking duplicates with samtools markdup:" >> {log} 2>&1
        echo "samtools markdup -r {output.sorted_bam} {output.bam_markdup}" >> {log} 2>&1
        (samtools markdup -r {output.sorted_bam} {output.bam_markdup}) >> {log} 2>&1

        echo "Indexing the sorted bam with marked duplicates file with samtools index" >> {log} 2>&1
        samtools index {output.bam_markdup} >> {log} 2>&1
        """  