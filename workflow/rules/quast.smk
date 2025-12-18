rule quast:
    input:
        assemblies=lambda wildcards: expand(
            f"{output_dir}assemblies/{{sample}}/{{sample}}_{{reads_type}}_{{strategy}}_{{assembler}}.fa",
            sample=wildcards.sample,  # Only this wildcard exists
            reads_type=READS_TYPES,   # Expand over all read types
            strategy=KMER_STRATEGIES, # Expand over all strategies
            assembler=ASSEMBLERS,     # Expand over all assemblers
        ),
    output:
        report=f"{output_dir}" + "quast/{sample}/report.txt"
    params:
        quast_dir=lambda wildcards, output: os.path.dirname(output.report),
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
        "logs/{sample}/quast.log",
    benchmark:
        "benchmark/{sample}/quast.txt"
    conda:
        "../envs/quast.yaml"
    container:
        "docker://quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2"
    shell:
        """
        echo "Running quast with the following command:" >> {log} 2>&1
        echo "quast {input.assemblies} -o {params.quast_dir} -t {threads} {params.optional_params}" >> {log} 2>&1
        quast {input.assemblies} -o {params.quast_dir} --min-contig 250 -t {threads} {params.optional_params} >> {log} 2>&1
        """
