rule select_best_assembly:
    input:
        busco_dirs=expand(
            f"{output_dir}" + "{{sample}}/busco_specific/{assembler}",
            assembler=ASSEMBLERS,
        ),
        quast_dir=f"{output_dir}" + "{sample}/quast",
    output:
        best_assemblies_dir=directory(f"{output_dir}" + "{sample}/best_assembly/"),
        assembly=f"{output_dir}" + "{sample}/best_assembly/{sample}_best_assembly.fa",
    params:
        sample="{sample}",
        results_dir=config["output_dir"],
    threads: 1
    log:
        "logs/{sample}/select_best_assembly.log",
    benchmark:
        "benchmark/{sample}/select_best_assembly.txt"
    conda:
        "../envs/basic.yaml"
    shell:
        """
        workflow/scripts/select_best_assembly.sh -s {params.sample} -r {params.results_dir} -o {output.best_assemblies_dir} >> {log} 2>&1
        """
