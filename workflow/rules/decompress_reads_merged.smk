# this rule decompresses gzipped FASTQ files for SparseAssembler
# it is used to prepare the input for SparseAssembler, which requires uncompressed FASTQ files
# the output is a temporary file that will be used in the next steps of the workflow
# the input is expected to be gzipped FASTQ files with a specific naming convention


rule decompress_reads:
    input:
        gz=f"{input_dir}" + "{sample}/{sample}_merge.fq.gz",
    output:
        fq=temp(f"{output_dir}" + "fqreads/{sample}/{sample}_merged.fq"),
    threads: 1
    resources:
        mem_mb=get_very_low_mem,
        partition=config["very_low"]["partition"],
    conda:
        "../envs/basic.yaml"
    log:
        "logs/{sample}/merged_decompress.log",
    benchmark:
        "benchmark/{sample}/merged_decompress.txt"
    shell:
        """
        echo "Decompressing {input.gz} to {output.fq}" >> {log} 2>&1
        zcat {input.gz} > {output.fq}
        echo "Done!" >> {log} 2>&1
        """
