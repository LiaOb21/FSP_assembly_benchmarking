rule busco_2:
    input:
        assembly=f"{output_dir}"
        + "assemblies/{sample}/{sample}_best_assembly_pypolca.fa",
        specific_db=f"{output_dir}" + "busco_db/{sample}/busco_db.txt",
    output:
        general_dir=directory(
            f"{output_dir}" + "best_assembly_qc/{sample}/busco_general_pypolca"
        ),
        specific_dir=directory(
            f"{output_dir}" + "best_assembly_qc/{sample}/busco_specific_pypolca"
        ),
    params:
        path_to_busco_bds=config["busco"]["path_to_busco_bds"],
        lineage_general=config["busco"]["lineage_general"],
        lineage_specific=lambda wildcards, input: open(input.specific_db).read().strip(),
        optional_params=" ".join(
            f"{k}" if v is True else f"{k} {v}"
            for k, v in config["busco"]["optional_params"].items()
            if v and v is not False and v != ""
        ),
    threads: get_medium_threads
    resources:
        mem_mb=get_medium_mem,
        partition=config["medium"]["partition"],
    log:
        "logs/{sample}/busco_best_assembly_pypolca.log",
    benchmark:
        "benchmark/{sample}/busco_best_assembly_pypolca.txt"
    conda:
        "../envs/busco.yaml"
    container:
        "docker://quay.io/biocontainers/busco:6.0.0--pyhdfd78af_1"
    shell:
        """
        echo "Running busco with the following command:" >> {log} 2>&1
        echo "busco -i {input.assembly} --out_path {output.general_dir} -l {params.path_to_busco_bds}{params.lineage_general} -f -c {threads} {params.optional_params}" >> {log} 2>&1
        busco -i {input.assembly} --out_path {output.general_dir} -l {params.path_to_busco_bds}{params.lineage_general} -f -c {threads} {params.optional_params} >> {log} 2>&1
        
        echo "Running busco with the following command:" >> {log} 2>&1
        echo "busco -i {input.assembly} --out_path {output.specific_dir} -l {params.path_to_busco_bds}{params.lineage_specific} -f -c {threads} {params.optional_params}" >> {log} 2>&1
        busco -i {input.assembly} --out_path {output.specific_dir} -l {params.path_to_busco_bds}{params.lineage_specific} -f -c {threads} {params.optional_params} >> {log} 2>&1
        """
