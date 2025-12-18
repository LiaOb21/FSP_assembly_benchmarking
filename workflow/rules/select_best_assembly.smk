rule select_best_assembly:
    input:
        busco_dirs=expand(
            f"{output_dir}" + "{strategy}/{{sample}}/busco_specific/{assembler}",
            assembler=ASSEMBLERS, strategy=KMER_STRATEGIES
        ),
        quast_dirs=expand(
            f"{output_dir}" + "{strategy}/{{sample}}/quast",
            strategy=KMER_STRATEGIES,  # All quast outputs
        ),
    output:
        best_assemblies_dir=directory(f"{output_dir}" + "best_assembly/{sample}/"),
        assembly=f"{output_dir}" + "best_assembly/{sample}/{sample}_best_assembly.fa",
    params:
        sample="{sample}",
        results_dir=output_dir,
    threads: 1
    resources:
        mem_mb=get_very_low_mem,
        partition=config["very_low"]["partition"],
    log:
        "logs/select_best_assembly/{sample}/select_best_assembly.log",
    benchmark:
        "benchmark/select_best_assembly/{sample}/select_best_assembly.txt"
    conda:
        "../envs/basic.yaml"
    container:
        "docker://debian:stable-slim"
    shell:
        """
        echo "Running the script as follows:" >> {log} 2>&1
        echo "workflow/scripts/select_best_assembly.sh -s {params.sample} -r {params.results_dir} -o {output.best_assemblies_dir}" >> {log} 2>&1
        workflow/scripts/select_best_assembly.sh -s {params.sample} -r {params.results_dir} -o {output.best_assemblies_dir} >> {log} 2>&1
        """
