# create a rule that aggregates all BUSCO results for a sample (this is because we need busco to be finished for all samples before we can select the best assembly)
rule aggregate_busco_results:
    input:
        busco_dirs=expand(f"{output_dir}" + "{{sample}}/busco_specific/{assembler}", assembler=ASSEMBLERS),
    output:
        marker=f"{output_dir}" + "{sample}/busco_complete.marker",
    threads: 1
    log:
        "logs/{sample}/aggregate_busco_results.log"  # Different log file
    conda:
        "../envs/basic.yaml"    
    shell:
        """
        touch {output.marker}
        """

rule select_best_assembly:
    input:
        busco_marker=f"{output_dir}" + "{sample}/busco_complete.marker",
        quast_dir=f"{output_dir}" + "{sample}/quast",
    output:
        best_assemblies_dir=directory(f"{output_dir}" + "best_assemblies/{sample}"),
        best_assembly_txt=f"{output_dir}" + "best_assemblies/{sample}/best_assembly.txt",
    params:
        sample="{sample}",
        results_dir=f"{output_dir}",
        output_dir=f"{output_dir}" + "best_assemblies",
    threads: 1
    log:
        "logs/{sample}/select_best_assembly.log"
    benchmark:
        "benchmark/{sample}/select_best_assembly.txt"
    conda:
        "../envs/basic.yaml"
    shell:
        """
        workflow/scripts/select_best_assembly.sh -s {params.sample} -r {params.results_dir} -o {params.output_dir} >> {log} 2>&1
        """