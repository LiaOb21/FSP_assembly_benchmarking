rule coverage_viz:
    input:
        coverage_stats=f"{output_dir}"
        + "{sample}/best_assembly_qc/bwa_mem2_samtools_pilon/{sample}_best_assembly_pilon_coverage_stats.txt",
        flagstat=f"{output_dir}"
        + "{sample}/best_assembly_qc/bwa_mem2_samtools_pilon/{sample}_best_assembly_pilon_flagstat.txt",
    output:
        coverage_plot=f"{output_dir}"
        + "{sample}/best_assembly_qc/coverage_viz_pilon/{sample}_best_assembly_pilon_coverage_plot.png",
        coverage_summary=f"{output_dir}"
        + "{sample}/best_assembly_qc/coverage_viz_pilon/{sample}_best_assembly_pilon_coverage_summary.txt",
    threads: 1
    resources:
        mem_mb=get_low_mem,
        partition=config["low"]["partition"],
    log:
        "logs/{sample}/coverage_viz_best_assembly_pilon.log",
    benchmark:
        "benchmark/{sample}/coverage_viz_best_assembly_pilon.txt"
    conda:
        "../envs/data_viz.yaml"
    shell:
        """
        echo "Creating coverage visualization for {wildcards.sample}" >> {log} 2>&1
        
        python workflow/scripts/coverage_visualisation.py \
            {input.coverage_stats} \
            {output.coverage_plot} \
            --sample {wildcards.sample} \
            --assembler "Best assembly improved with Pilon" \
            --flagstat {input.flagstat} \
            --summary {output.coverage_summary} \
            >> {log} 2>&1
        
        echo "Coverage visualization completed successfully" >> {log} 2>&1
        """
