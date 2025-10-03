rule busco:
    input:
        assembly=f"{output_dir}" + "{sample}/assemblies/{sample}_{assembler}.fa",
        specific_db=f"{output_dir}" + "{sample}/busco_specific/busco_db.txt",
    output:
        general_dir=directory(f"{output_dir}" + "{sample}/busco_general/{assembler}"),
        specific_dir=directory(f"{output_dir}" + "{sample}/busco_specific/{assembler}"),
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
        "logs/{sample}/busco_{sample}_{assembler}.log",
    benchmark:
        "benchmark/{sample}/busco_{sample}_{assembler}.txt"
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
