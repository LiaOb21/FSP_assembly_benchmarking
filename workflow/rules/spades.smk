# This rule runs spades assembler with the parameters specified in the config file.
import glob


rule spades:
    input:
        forward_in= "{config['input_path']}/{sample}/{sample}_trimmed.R1.fq.gz",
        reverse_in= "{config['input_path']}/{sample}/{sample}_trimmed.R2.fq.gz",
    output:
        result_dir = directory("{config['results_path']}/spades/{sample}"),
        scaffolds = "{config['results_path']}/spades/{sample}/scaffolds.fasta"
    params:
        optional_params=" ".join(
            f"{k} {v}" for k, v in config["spades"]["optional_params"].items() if v
        ),
    threads: config["threads"],  # access threads from config
    log:
        "logs/{sample}/spades.log",
    resources:
        mem_mb=config["mem_mb"],  # access memory from config
    conda:
        "../envs/spades.yaml"
    message:
        "Running: spades.py -t {threads} -1 {input.forward_in} -2 {input.reverse_in} -o {output.result_dir} {params.optional_params} >> {log} 2>&1"
    
    run:
        print(f"[spades] Input: {input.forward_in}, {input.reverse_in}")
        print(f"[spades] Output: {output.scaffolds}")
    
    #shell:
    #    """
    #    spades.py -t {threads} -1 {input.forward_in} -2 {input.reverse_in} -o {output.result_dir} {params.optional_params} >> {log} 2>&1
    #    """