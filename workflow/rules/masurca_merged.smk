# This rule runs masurca assembler with the parameters specified in the config file.


rule masurca_merged:
    input:
        merged_in=f"{input_dir}" + "{sample}/{sample}_merge.fq.gz",
        template="workflow/scripts/masurca_config_template_merged.txt",
        kmergenie_result=get_kmergenie_dependency,
    output:
        masurca_config=f"{output_dir}" + "{sample}/masurca/masurca_config.txt",
        scaffolds=f"{output_dir}" + "{sample}/masurca/CA/primary.genome.scf.fasta",
        link_assembly=f"{output_dir}" + "{sample}/assemblies/{sample}_masurca.fa",
    params:
        fragment_mean=config["masurca"].get("fragment_mean", 500),
        fragment_stdev=config["masurca"].get("fragment_stdev", 50),
        k=lambda wildcards: get_single_kmer(wildcards, "masurca", "k"),
        jf_size=config["masurca"].get("jf_size", 10000000000),
        ca_parameters=config["masurca"].get("ca_parameters", "cgwErrorRate=0.15"),
        config_dir=lambda wildcards, input, output: os.path.dirname(
            output.masurca_config
        ),
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
    container:
        "docker://quay.io/biocontainers/masurca:4.1.4--h6b3f7d6_0"
    shell:
        """
        echo "Processing sample: {wildcards.sample}" >> {log} 2>&1
        echo "Using k-mer size: {params.k}" >> {log} 2>&1 
        python3 <<EOF
with open("{input.template}") as t, open("{output.masurca_config}", "w") as out:
    template = t.read()
    out.write(template.format(
        input_merged="{input.merged_in}",
        fragment_mean={params.fragment_mean},
        fragment_stdev={params.fragment_stdev},
        kmer="{params.k}",
        threads={threads},
        jf_size={params.jf_size},
        ca_parameters="{params.ca_parameters}"
    ))
EOF
        echo "Config generated for sample {wildcards.sample} with kmer size {params.k}" >> {log} 2>&1

        cd {params.config_dir}
        masurca {output.masurca_config} 
        ./assemble.sh >> masurca.log 2>&1

        cd -

        cp {params.config_dir}/masurca.log {log}
        
        ln -srn {output.scaffolds} {output.link_assembly}
        """


# Leaving as a note:
# masurca results/048ss/masurca/masurca_config.txt -o results/048ss/masurca/assemble.sh
# sed -e "s/memory 1000000000/memory 1000000/g" results/048ss/masurca/assemble.sh > results/048ss/masurca/assemble_1.sh
# ./results/048ss/masurca/assemble_1.sh this unfortunately writes in the current directory. Also it's still failing with core dumped
