rule masurca_config:
    input:
        merged_in=f"{input_dir}" + "{sample}/{sample}_merge.fq.gz",
        template="workflow/scripts/masurca_config_template_merged.txt",
        kmergenie_result=get_kmergenie_dependency,
    output:
        cfg=protected(f"{output_dir}" + "{sample}/masurca/masurca_config.txt"),
    threads: get_high_threads
    resources:
        mem_mb=get_low_mem,
    params:
        fragment_mean=config["masurca"].get("fragment_mean", 500),
        fragment_stdev=config["masurca"].get("fragment_stdev", 50),
        k=lambda wildcards: get_single_kmer(wildcards, "masurca", "k"),
        jf_size=config["masurca"].get("jf_size", 10000000000),
        ca_parameters=config["masurca"].get("ca_parameters", "cgwErrorRate=0.15"),
    log:
        "logs/{sample}/masurca_config.log",
    benchmark:
        "benchmark/{sample}/masurca_config.txt"
    conda:
        "../envs/basic.yaml"
    shell:
        """
        echo "Processing sample: {wildcards.sample}" >> {log} 2>&1
        echo "Using k-mer size: {params.k}" >> {log} 2>&1 
        python3 <<EOF
with open("{input.template}") as t, open("{output.cfg}", "w") as out:
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
        """
