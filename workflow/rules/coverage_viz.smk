rule coverage_viz:
    input:
        coverage_stats=f"{output_dir}"
        + "best_assembly_qc/{sample}/samtools_pypolca/{sample}_best_assembly_pypolca_coverage_stats.txt",
        flagstat=f"{output_dir}"
        + "best_assembly_qc/{sample}/samtools_pypolca/{sample}_best_assembly_pypolca_flagstat.txt",
    output:
        coverage_plot=f"{output_dir}"
        + "best_assembly_qc/{sample}/coverage_viz_pypolca/{sample}_best_assembly_pypolca_coverage_plot.png",
        coverage_summary=f"{output_dir}"
        + "best_assembly_qc/{sample}/coverage_viz_pypolca/{sample}_best_assembly_pypolca_coverage_summary.txt",
    threads: 1
    resources:
        mem_mb=get_very_low_mem,
        partition=config["very_low"]["partition"],
    log:
        "logs/coverage_viz_best_assembly/{sample}/coverage_viz_best_assembly_pypolca.log",
    benchmark:
        "benchmark/{sample}/coverage_viz_best_assembly_pypolca.txt"
    conda:
        "../envs/data_viz.yaml"
    container:
        "docker://lobinu21/data_viz_py:latest"
    shell:
        """
        echo "Creating coverage visualization for {wildcards.sample}" >> {log} 2>&1
        
        python workflow/scripts/coverage_visualisation.py \
            {input.coverage_stats} \
            {output.coverage_plot} \
            --sample {wildcards.sample} \
            --assembler "Best assembly improved with pypolca" \
            --flagstat {input.flagstat} \
            --summary {output.coverage_summary} \
            >> {log} 2>&1
        
        echo "Coverage visualization completed successfully" >> {log} 2>&1
        """
