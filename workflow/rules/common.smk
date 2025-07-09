# resource scaling functions


def get_scaled_mem(wildcards, attempt):
    base_mem = config["mem_mb"]
    return min(base_mem * attempt, 250000)  # Cap at 250GB


def get_scaled_threads(wildcards, attempt):
    base_threads = config["threads"]
    return min(
        base_threads + (attempt - 1) * 4, 64
    )  # Add 4 threads per retry, cap at 64
