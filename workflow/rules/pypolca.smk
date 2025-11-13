rule pypolca:
    input:
        forward_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R1.fq.gz",
        reverse_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R2.fq.gz",
        assembly=f"{output_dir}" + "{sample}/best_assembly/{sample}_best_assembly.fa",
    output:
        pypolca_fasta=f"{output_dir}"
        + "{sample}/best_assembly/pypolca/{sample}_best_assembly_pypolca_corrected.fasta",
        link_pypolca_assembly=f"{output_dir}"
        + "{sample}/assemblies/{sample}_best_assembly_pypolca.fa",
    params:
        results_prefix=lambda wildcards, output: os.path.splitext(output.pypolca_fasta)[
            0
        ].replace(".fasta", ""),
        out_dir=lambda wildcards, input, output: os.path.dirname(output.pypolca_fasta),
        optional_params=" ".join(
            f"{key}" if value is True else f"{key} {value}"
            for key, value in config["pypolca"]["optional_params"].items()
            if value and value is not False and value != ""
        ),
    threads: get_high_threads
    resources:
        mem_mb=get_high_mem,
        partition=config["high"]["partition"],
    log:
        "logs/{sample}/pypolca_best_assembly.log",
    benchmark:
        "benchmark/{sample}/pypolca_best_assembly.txt"
    conda:
        "../envs/pypolca.yaml"
    shell:
        """
        echo "Running pypolca with the following command:" >> {log} 2>&1
        echo "pypolca run -a {input.assembly} -1 {input.forward_in} -2 {input.reverse_in} -t {threads} -f -o {params.out_dir} -p {params.results_prefix} {params.optional_params}" >> {log} 2>&1
        pypolca run -a {input.assembly} -1 {input.forward_in} -2 {input.reverse_in} -t {threads} -f -o {params.out_dir} -p {params.results_prefix} {params.optional_params} >> {log} 2>&
        
        ln -srn {output.pypolca_fasta} {output.link_pypolca_assembly}
        """
