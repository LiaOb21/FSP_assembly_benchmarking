# This rule runs pypolca on the best draft assemblies


rule pypolca:
    input:
        forward_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R1.fq.gz",
        reverse_in=f"{input_dir}" + "{sample}/{sample}_trimmed.R2.fq.gz",
        assembly=f"{output_dir}" + "best_assembly/{sample}/{sample}_best_assembly.fa",
    output:
        pypolca_fasta=f"{output_dir}"
        + "best_assembly/{sample}/pypolca/{sample}_best_assembly_pypolca_corrected.fasta",
        link_pypolca_assembly=f"{output_dir}"
        + "best_assembly_fa/{sample}/{sample}_best_assembly_pypolca.fa",
    params:
        results_prefix=lambda wildcards, output: os.path.basename(
            output.pypolca_fasta
        ).replace("_corrected.fasta", ""),
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
        "logs/pypolca/{sample}/pypolca_best_assembly.log",
    benchmark:
        "benchmark/pypolca/{sample}/pypolca_best_assembly.txt"
    conda:
        "../envs/pypolca.yaml"
    container:
        "docker://quay.io/biocontainers/pypolca:0.4.0--pyhdfd78af_0"
    shell:
        """
        echo "Running pypolca with the following command:" >> {log} 2>&1
        echo "pypolca run -a {input.assembly} -1 {input.forward_in} -2 {input.reverse_in} -t {threads} -f -o {params.out_dir} -p {params.results_prefix} {params.optional_params}" >> {log} 2>&1
        pypolca run -a {input.assembly} -1 {input.forward_in} -2 {input.reverse_in} -t {threads} -f -o {params.out_dir} -p {params.results_prefix} {params.optional_params} >> {log} 2>&1
        
        ln -srn {output.pypolca_fasta} {output.link_pypolca_assembly}
        """
