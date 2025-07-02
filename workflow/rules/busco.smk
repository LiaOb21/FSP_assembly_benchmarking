rule busco:
    input:
        assembly=f"{output_dir}" + "assemblies/{sample}/{sample}_{assembler}.fa",
    output:
        dir=directory(f"{output_dir}" + "{sample}/busco/{assembler}"),
    params:
        lineage=config["busco"]["lineage"],
        optional_params=" ".join(
            f"{k}" if v is True else f"{k} {v}"
            for k, v in config["busco"]["optional_params"].items()
            if v
        ),
    threads: config["threads"]
    log:
        "logs/{sample}/busco_{sample}_{assembler}.log",
    resources:
        mem_mb=config["mem_mb"],
    conda:
        "../envs/busco.yaml"
    shell:
        """
        busco -i {input.assembly} --out_path {output.dir} -l {params.lineage} -f -c {threads} {params.optional_params} >> {log} 2>&1
        """
