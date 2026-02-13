rule select_best_assembly:
    input:
        busco_dirs=expand(
            f"{output_dir}"
            + "{reads_type}/{strategy}/{{sample}}/busco_specific/{assembler}",
            reads_type=READS_TYPES,
            strategy=KMER_STRATEGIES,
            assembler=ASSEMBLERS,
        ),
        report=f"{output_dir}" + "quast/{sample}/report.txt",
    output:
        assembly=f"{output_dir}" + "best_assembly/{sample}/{sample}_best_assembly.fa",
    params:
        sample="{sample}",
        results_dir=lambda w, output: os.path.dirname(
            os.path.dirname(os.path.dirname(output.assembly))
        )
        + "/",
        best_assemblies_dir=lambda w, output: os.path.dirname(output.assembly)
        + "/",
    threads: 1
    resources:
        mem_mb=get_very_low_mem,
        partition=config["very_low"]["partition"],
    log:
        "logs/{sample}/select_best_assembly.log",
    benchmark:
        "benchmark/{sample}/select_best_assembly.txt"
    conda:
        "../envs/basic.yaml"
    container:
        "docker://debian:stable-slim"
    shell:
        """
        echo "Running the script as follows:" >> {log} 2>&1
        echo "workflow/scripts/select_best_assembly.sh -s {params.sample} -r {params.results_dir} -o {params.best_assemblies_dir}" >> {log} 2>&1
        workflow/scripts/select_best_assembly.sh -s {params.sample} -r {params.results_dir} -o {params.best_assemblies_dir} >> {log} 2>&1
        """
