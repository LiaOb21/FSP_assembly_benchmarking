rule quast:
    input:
        assemblies=lambda wildcards: [
            f"{output_dir}assemblies/{wildcards.sample}/{wildcards.sample}_{assembler}.fa"
            for assembler in ASSEMBLERS
        ]
    output:
        dir = directory(f"{output_dir}" + "{sample}/quast"),
    params:
        optional_params = " ".join(
            f"{k}" if v is True else f"{k} {v}"
            for k, v in config["quast"]["optional_params"].items() if v
        ),
    threads: config["threads"]
    log:
        "logs/{sample}/quast.log"
    resources:
        mem_mb = config["mem_mb"]
    conda:
        "../envs/quast.yaml"
    shell:
        """
        quast {input.assemblies} -o {output.dir} -t {threads} {params.optional_params} >> {log} 2>&1
        """