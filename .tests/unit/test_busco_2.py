import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_busco_2():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/busco_2/data")
        expected_path = PurePosixPath(".tests/unit/busco_2/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print("/home/lobinu/git_repos/FSP_assembly_benchmarking/results_3/best_assembly_qc/048ds/busco_general_pypolca /home/lobinu/git_repos/FSP_assembly_benchmarking/results_3/best_assembly_qc/048ds/busco_specific_pypolca", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "/home/lobinu/git_repos/FSP_assembly_benchmarking/results_3/best_assembly_qc/048ds/busco_general_pypolca /home/lobinu/git_repos/FSP_assembly_benchmarking/results_3/best_assembly_qc/048ds/busco_specific_pypolca",
            "-f", 
            "-j1",
            "--target-files-omit-workdir-adjustment",
            "--configfile",
            /home/lobinu/git_repos/FSP_assembly_benchmarking/config/test_config.yml
    
            "--directory",
            workdir,
        ])

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
