# this rule decompresses gzipped FASTQ files for SparseAssembler
# it is used to prepare the input for SparseAssembler, which requires uncompressed FASTQ files
# the output is a temporary file that will be used in the next steps of the workflow
# the input is expected to be gzipped FASTQ files with a specific naming convention

rule decompress_reads:
    input:
        gz = f"{input_dir}" + "/{sample}/{sample}_trimmed.{read}.fq.gz"
    output:
        fq = temp(f"{output_dir}" + "/fqreads/{sample}/{sample}_trimmed.{read}.fq")
    shell:
        """
        zcat {input.gz} > {output.fq}
        """