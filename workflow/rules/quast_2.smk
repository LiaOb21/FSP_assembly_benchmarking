rule quast_2:
    input:
        assembly=f"{output_dir}" + "{sample}/assemblies/{sample}_best_assembly_pilon.fa",
    output:
        dir=directory(f"{output_dir}" + "{sample}/best_assembly_qc/quast_pilon"),
    params:
        optional_params=" ".join(
            f"{k}" if v is True else f"{k} {v}"
            for k, v in config["quast"]["optional_params"].items()
            if v and v is not False and v != ""
        ),
    threads: get_scaled_threads  # Use scaling function
    log:
        "logs/{sample}/quast_best_assembly_pilon.log",
    benchmark:
        "benchmark/{sample}/quast_best_assembly_pilon.txt"
    resources:
        mem_mb=get_scaled_mem,  # Use scaling function
    conda:
        "../envs/quast.yaml"
    shell:
        """
        echo "Running quast with the following command:" >> {log} 2>&1
        echo "quast {input.assembly} -o {output.dir} -t {threads} {params.optional_params}" >> {log} 2>&1
        quast {input.assembly} -o {output.dir} --min-contig 250 -t {threads} {params.optional_params} >> {log} 2>&1
        """
