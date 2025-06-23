# This rule runs masurca assembler with the parameters specified in the config file.
import glob
import os


rule masurca:
    input:
        masurca_config= f"{output_dir}" + "/{sample}/masurca/masurca_config.txt",
    output:
        assemble_sh = directory(f"{output_dir}" + "/{sample}/masurca"),
        scaffolds = f"{output_dir}" + "/{sample}/masurca/CA/primary.genome.scf.fasta",
    threads: config["threads"],  # access threads from config
    params:
#        config_dir = f"{output_dir}" + "/{sample}/masurca",
    log:
        "logs/{sample}/masurca.log",
    resources:
        mem_mb=config["mem_mb"],  # access memory from config
    conda:
        "../envs/masurca.yaml"
    shell:
        """
        cd {params.config_dir}
        masurca {input.masurca_config} >> masurca_run.log 2>&1
        mv assemble.sh {config_dir}
        sed -i -E "s/memory 1000000000/memory {resources.mem_mb}00000/g" assemble.sh > assemble_1.sh
        ./assemble_1.sh >> masurca_run.log 2>&1

        """

#Leaving as a note: 

#masurca results/048ss/masurca/masurca_config.txt -o results/048ss/masurca/assemble.sh
#sed -e "s/memory 1000000000/memory 1000000/g" results/048ss/masurca/assemble.sh > results/048ss/masurca/assemble_1.sh
# ./results/048ss/masurca/assemble_1.sh this unfortunately writes in the current directory. Also it's still failing with core dumped