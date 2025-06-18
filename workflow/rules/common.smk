
def get_all_inputs(wc=None):
    return [
        f"{config['results_path']}/spades/{sample}/scaffolds.fasta"
        for sample in samples
    ]
