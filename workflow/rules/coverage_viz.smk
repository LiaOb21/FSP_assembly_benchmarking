rule coverage_viz:
    input:
        coverage_stats=f"{output_dir}" + "{sample}/bwa_mem2_samtools/{assembler}/{sample}_{assembler}_coverage_stats.txt",
        flagstat=f"{output_dir}" + "{sample}/bwa_mem2_samtools/{assembler}/{sample}_{assembler}_flagstat.txt",
    output:
        coverage_plot=f"{output_dir}" + "{sample}/coverage_viz/{assembler}/{sample}_{assembler}_coverage_plot.png",
        coverage_summary=f"{output_dir}" + "{sample}/coverage_viz/{assembler}/{sample}_{assembler}_coverage_summary.txt",
    threads: 1
    log:
        "logs/{sample}/coverage_viz_{assembler}.log",
    benchmark:
        "benchmark/{sample}/coverage_viz_{assembler}.txt"
    conda:
        "../envs/data_viz.yaml"
    shell:
        """
        echo "Creating coverage visualization for {wildcards.sample} - {wildcards.assembler}" >> {log} 2>&1
        
        python workflow/scripts/coverage_visualisation.py \
            {input.coverage_stats} \
            {output.coverage_plot} \
            --sample {wildcards.sample} \
            --assembler {wildcards.assembler} \
            --flagstat {input.flagstat} \
            --summary {output.coverage_summary} \
            >> {log} 2>&1
        
        echo "Coverage visualization completed successfully" >> {log} 2>&1
        """