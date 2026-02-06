# this rule decompresses gzipped FASTQ files for SparseAssembler
# it is used to prepare the input for SparseAssembler, which requires uncompressed FASTQ files
# the output is a temporary file that will be used in the next steps of the workflow
# the input is expected to be gzipped FASTQ files with a specific naming convention


rule decompress_reads:
    wildcard_constraints:
        reads_type="R1R2",
    input:
        gz=f"{input_dir}" + "{sample}/{sample}_trimmed.{read}.fq.gz",
    output:
        fq=temp(
            f"{output_dir}"
            + "{reads_type}/fqreads/{sample}/{sample}_trimmed.{read}.fq"
        ),
    threads: 1
    resources:
        mem_mb=get_very_low_mem,
        partition=config["very_low"]["partition"],
    conda:
        "../envs/basic.yaml"
    container:
        "docker://debian:stable-slim"
    log:
        "logs/{sample}/{reads_type}_{read}_decompress.log",
    benchmark:
        "benchmark/{sample}/{reads_type}_{read}_decompress.txt"
    shell:
        """
        echo "Decompressing {input.gz} to {output.fq}" >> {log} 2>&1
        zcat {input.gz} > {output.fq}
        echo "Done!" >> {log} 2>&1
        """
