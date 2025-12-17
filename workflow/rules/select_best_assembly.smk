rule select_best_assembly:
    input:
        busco_dirs=expand(
            f"{output_dir}" + "{strategy}/{{sample}}/busco_specific/{assembler}",
            assembler=ASSEMBLERS, strategy=KMER_STRATEGIES
        ),
        quast_dir=f"{output_dir}" + "{strategy}/{sample}/quast",
    output:
        best_assemblies_dir=directory(f"{output_dir}" + "{strategy}/{sample}/best_assembly/"),
        assembly=f"{output_dir}" + "{strategy}/{sample}/best_assembly/{sample}_best_assembly.fa",
    params:
        sample="{sample}",
        results_dir=lambda wildcards, input: os.path.dirname(
            os.path.dirname(input.quast_dir)
        ),
    threads: 1
    resources:
        mem_mb=get_very_low_mem,
        partition=config["very_low"]["partition"],
    log:
        "logs/{strategy}/{sample}/select_best_assembly.log",
    benchmark:
        "benchmark/{strategy}/{sample}/select_best_assembly.txt"
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
