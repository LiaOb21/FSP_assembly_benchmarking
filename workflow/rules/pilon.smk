# This rule runs pilon on the best draft assemblies
import glob
import os


rule pilon:
    input:
        assembly=f"{output_dir}" + "{sample}/best_assembly/{sample}_best_assembly.fa",
        bam_markdup=f"{output_dir}"
        + "{sample}/best_assembly/samtools/{sample}_best_assembly_fixmate_sorted_markdup.bam",
    output:
        pilon_fasta=f"{output_dir}"
        + "{sample}/best_assembly/pilon/{sample}_best_assembly_pilon.fasta",
        link_pilon_assembly=f"{output_dir}"
        + "{sample}/assemblies/{sample}_best_assembly_pilon.fa",
    params:
        results_prefix=lambda wildcards, output: os.path.splitext(output.pilon_fasta)[
            0
        ].replace(".fasta", ""),
        optional_params=" ".join(
            f"{key}" if value is True else f"{key} {value}"
            for key, value in config["pilon"]["optional_params"].items()
            if value and value is not False and value != ""
        ),
        java_heap=lambda wildcards, resources: f"{int(resources.mem_mb*0.8//1024)}G",
    threads: get_medium_threads
    resources:
        mem_mb=get_high_mem,
        partition=config["high"]["partition"],
    log:
        "logs/{sample}/pilon_best_assembly.log",
    benchmark:
        "benchmark/{sample}/pilon_best_assembly.txt"
    conda:
        "../envs/pilon.yaml"
    shell:
        """
        PILON=${{CONDA_PREFIX}}/share/pilon-*/pilon*.jar
        echo "Running pilon with the following command:" >> {log} 2>&1
        echo "java -Xmx{params.java_heap} -jar ${{PILON}} --genome {input.assembly} --bam {input.bam_markdup} --output {params.results_prefix} --threads {threads} {params.optional_params} --verbose" >> {log} 2>&1
        java -Xmx{params.java_heap} -jar ${{PILON}} --genome {input.assembly} --bam {input.bam_markdup} --output {params.results_prefix} --threads {threads} {params.optional_params} --verbose >> {log} 2>&1
        ln -srn {output.pilon_fasta} {output.link_pilon_assembly}
        """
