# This rule runs masurca assembler with the parameters specified in the config file.
import glob
import os


rule masurca:
    input:
        masurca_config=f"{output_dir}" + "{sample}/masurca/masurca_config.txt",
    output:
        scaffolds=f"{output_dir}" + "{sample}/masurca/CA/primary.genome.scf.fasta",
        link_assembly=f"{output_dir}" + "{sample}/assemblies/{sample}_masurca.fa",
    params:
        config_dir=lambda wildcards, input: os.path.dirname(input.masurca_config),
    threads: get_high_threads
    resources:
        mem_mb=get_high_mem,
        partition=config["high"]["partition"],
    log:
        "logs/{sample}/masurca.log",
    benchmark:
        "benchmark/{sample}/masurca.txt"
    conda:
        "../envs/masurca.yaml"
    shell:
        """
        cd {params.config_dir}
        masurca {input.masurca_config} 
        ./assemble.sh >> masurca.log 2>&1

        cd -

        cp {params.config_dir}/masurca.log {log}
        
        ln -srn {output.scaffolds} {output.link_assembly}
        """


# Leaving as a note:
# masurca results/048ss/masurca/masurca_config.txt -o results/048ss/masurca/assemble.sh
# sed -e "s/memory 1000000000/memory 1000000/g" results/048ss/masurca/assemble.sh > results/048ss/masurca/assemble_1.sh
# ./results/048ss/masurca/assemble_1.sh this unfortunately writes in the current directory. Also it's still failing with core dumped
