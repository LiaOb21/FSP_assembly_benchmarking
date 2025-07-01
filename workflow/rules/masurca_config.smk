rule masurca_config:
    input:
        r1 = f"{input_dir}" + "/{sample}/{sample}_trimmed.R1.fq.gz",
        r2 = f"{input_dir}" + "/{sample}/{sample}_trimmed.R2.fq.gz",
        template = "workflow/scripts/masurca_config_template.txt"
    output:
        cfg = protected(f"{output_dir}" + "/{sample}/masurca/masurca_config.txt")
    params:
        fragment_mean = config["masurca"].get("fragment_mean", 500),
        fragment_stdev = config["masurca"].get("fragment_stdev", 50),
        kmer = config["masurca"].get("kmer", "auto"),
        threads = config["threads"],
        jf_size = config["masurca"].get("jf_size", 10000000000),
        ca_parameters = config["masurca"].get("ca_parameters", "cgwErrorRate=0.15"),
    shell:
        """
        python3 <<EOF
with open("{input.template}") as t, open("{output.cfg}", "w") as out:
    template = t.read()
    out.write(template.format(
        input_r1="{input.r1}",
        input_r2="{input.r2}",
        fragment_mean={params.fragment_mean},
        fragment_stdev={params.fragment_stdev},
        kmer="{params.kmer}",
        threads={params.threads},
        jf_size={params.jf_size},
        ca_parameters="{params.ca_parameters}"
    ))
EOF
        """