rule quast_2:
    input:
        assembly=f"{output_dir}"
        + "{strategy}/{sample}/assemblies/{sample}_best_assembly_pypolca.fa",
    output:
        dir=directory(f"{output_dir}" + "{strategy}/{sample}/best_assembly_qc/quast_pypolca"),
    params:
        optional_params=" ".join(
            f"{k}" if v is True else f"{k} {v}"
            for k, v in config["quast"]["optional_params"].items()
            if v and v is not False and v != ""
        ),
    threads: get_low_threads
    resources:
        mem_mb=get_low_mem,
        partition=config["low"]["partition"],
    log:
        "logs/{strategy}/{sample}/quast_best_assembly_pypolca.log",
    benchmark:
        "benchmark/{strategy}/{sample}/quast_best_assembly_pypolca.txt"
    conda:
        "../envs/quast.yaml"
    container:
        "docker://quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2"
    shell:
        """
        echo "Running quast with the following command:" >> {log} 2>&1
        echo "quast {input.assembly} -o {output.dir} -t {threads} {params.optional_params}" >> {log} 2>&1
        quast {input.assembly} -o {output.dir} --min-contig 250 -t {threads} {params.optional_params} >> {log} 2>&1
        """
