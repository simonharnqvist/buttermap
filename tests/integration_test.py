import pytest
import secrets
import shutil
import subprocess
from pathlib import Path

def test_run_pipeline(tmp_path):

    print(tmp_path)
    subprocess.run([f"cp test_data/test_input/* {tmp_path}"], shell=True, check=True)

    subprocess.run(["buttermap", 
                "--region-size", "100000",
                "--reference-fasta", f"{tmp_path}/brenthis_ino_chr14.reference.fasta", 
                "--reads-dir", f"{tmp_path}",
                "--output-dir", f"{tmp_path}"], check=True)
    
    correct_vcf_name = "brenthis_ino_chr14.vcf.gz"
    assert Path(f"{tmp_path}/{correct_vcf_name}").exists()

    subprocess.run(["bgzip", "-d", f"{tmp_path}/{correct_vcf_name}"], check=True)
    with open(f"{tmp_path}/brenthis_ino_chr14.vcf", "r") as f:
        line0 = f.readlines()[0]
    assert line0 == "##fileformat=VCFv4.2\n"