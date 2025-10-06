rule busco_2:
    input:
        assembly=f"{output_dir}" + "{sample}/assemblies/{sample}_best_assembly_pilon.fa",
        specific_db=f"{output_dir}" + "{sample}/busco_db/busco_db.txt",
    output:
        general_dir=directory(
            f"{output_dir}" + "{sample}/best_assembly_qc/busco_general_pilon"
        ),
        specific_dir=directory(
            f"{output_dir}" + "{sample}/best_assembly_qc/busco_specific_pilon"
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
        "logs/{sample}/busco_best_assembly_pilon.log",
    benchmark:
        "benchmark/{sample}/busco_best_assembly_pilon.txt"
    conda:
        "../envs/busco.yaml"
    shell:
        """
        echo "Running busco with the following command:" >> {log} 2>&1
        echo "busco -i {input.assembly} --out_path {output.general_dir} -l {params.path_to_busco_bds}{params.lineage_general} -f -c {threads} {params.optional_params}" >> {log} 2>&1
        busco -i {input.assembly} --out_path {output.general_dir} -l {params.path_to_busco_bds}{params.lineage_general} -f -c {threads} {params.optional_params} >> {log} 2>&1
        
        echo "Running busco with the following command:" >> {log} 2>&1
        echo "busco -i {input.assembly} --out_path {output.specific_dir} -l {params.path_to_busco_bds}{params.lineage_specific} -f -c {threads} {params.optional_params}" >> {log} 2>&1
        busco -i {input.assembly} --out_path {output.specific_dir} -l {params.path_to_busco_bds}{params.lineage_specific} -f -c {threads} {params.optional_params} >> {log} 2>&1
        """
