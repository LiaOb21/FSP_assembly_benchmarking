# This rule generates a runs merquryFK to evaluate assembly quality based on kmer


rule merquryfk_2:
    input:
        ktab=get_fastk_table_for_best_assembly,
        assembly=f"{output_dir}"
        + "assemblies/{sample}/{sample}_best_assembly_pypolca.fa",
    output:
        stats=f"{output_dir}"
        + "best_assembly_qc/{sample}/merquryfk_pypolca/merquryfk.completeness.stats",
        qv=f"{output_dir}" + "best_assembly_qc/{sample}/merquryfk_pypolca/merquryfk.qv",
    params:
        result_prefix=lambda wildcards, output: os.path.splitext(output.qv)[0],
        temp_dir=lambda wildcards, output: os.path.join(
            os.path.dirname(output.qv), "temp"
        ),
        optional_params=" ".join(
            f"{k}" if v is True else f"{k} {v}"
            for k, v in config["merquryfk"]["optional_params"].items()
            if v and v is not False and v != ""
        ),
    threads: get_low_threads
    resources:
        mem_mb=get_low_mem,
        partition=config["low"]["partition"],
    log:
        "logs/{sample}/merquryfk_best_assembly_pypolca.log",
    benchmark:
        "benchmark/{sample}/merquryFK_best_assembly_pypolca.txt"
    conda:
        "../envs/merquryFK.yaml"
    container:
        "docker://quay.io/biocontainers/merquryfk:1.1.3--h71df26d_0"
    shell:
        """
        mkdir -p {params.temp_dir}

        echo "Running merquryFK with the following command:" >> {log} 2>&1
        echo "MerquryFK -lfs -v -T{threads} -P{params.temp_dir} {input.ktab} {input.assembly} {params.result_prefix} {params.optional_params}" >> {log} 2>&1
        MerquryFK -lfs -v -T{threads} -P{params.temp_dir} {input.ktab} {input.assembly} {params.result_prefix} {params.optional_params} >> {log} 2>&1
        rm -rf {params.temp_dir}
        """
