# This rule indexes the assemblies produced by the workflow, alignes the raw reads to the genome, and produces statistics


rule bwa_mem2:
    input:
        assembly=f"{output_dir}"
        + "assemblies/{sample}/{sample}_best_assembly_pypolca.fa",
        forward_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R1.fq.gz",
        reverse_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R2.fq.gz",
    output:
        sam=f"{output_dir}"
        + "best_assembly_qc/{sample}/bwa_mem2_pypolca/{sample}_best_assembly_pypolca.sam",
    threads: get_medium_high_threads
    resources:
        mem_mb=get_medium_high_mem,
        partition=config["medium_high"]["partition"],
    log:
        "logs/{sample}/bwa_mem2_best_assembly_pypolca.log",
    benchmark:
        "benchmark/{sample}/bwa_mem2_best_assembly_pypolca.txt"
    conda:
        "../envs/bwa_mem2.yaml"
    container:
        "docker://quay.io/biocontainers/bwa-mem2:2.3--he70b90d_0"
    shell:
        """
        echo "Running bwa-mem2 index with the following command:" >> {log} 2>&1
        echo "bwa-mem2 index {input.assembly}" >> {log} 2>&1
        bwa-mem2 index {input.assembly} >> {log} 2>&1

        echo "Running bwa-mem2 mem with the following command:" >> {log} 2>&1
        echo "bwa-mem2 mem -t {threads} {input.assembly} {input.forward_in} {input.reverse_in} > {output.sam}" >> {log} 2>&1
        (bwa-mem2 mem -t {threads} {input.assembly} {input.forward_in} {input.reverse_in} > {output.sam}) >> {log} 2>&1
        """
