rule quast:
    input:
        assemblies=lambda wildcards: [
            f"{output_dir}assemblies/{wildcards.sample}/{wildcards.sample}_{assembler}.fa"
            for assembler in ASSEMBLERS
        ],
        assemblies_pilon=lambda wildcards: [
            f"{output_dir}assemblies/{wildcards.sample}/{wildcards.sample}_{assembler}_pilon.fa"
            for assembler in ASSEMBLERS
        ],
    output:
        dir=directory(f"{output_dir}" + "{sample}/quast"),
    params:
        optional_params=" ".join(
            f"{k}" if v is True else f"{k} {v}"
            for k, v in config["quast"]["optional_params"].items()
            if v and v is not False and v != ""
        ),
    threads: get_scaled_threads  # Use scaling function
    log:
        "logs/{sample}/quast.log",
    benchmark:
        "benchmark/{sample}/quast.txt"
    resources:
        mem_mb=get_scaled_mem,  # Use scaling function
    conda:
        "../envs/quast.yaml"
    shell:
        """
        echo "Running quast with the following command:" >> {log} 2>&1
        echo "quast {input.assemblies} {input.assemblies_pilon} -o {output.dir} -t {threads} {params.optional_params}" >> {log} 2>&1
        quast {input.assemblies} {input.assemblies_pilon} -o {output.dir} -t {threads} {params.optional_params} >> {log} 2>&1
        """
