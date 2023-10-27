from buttermap.freebayes_region import run_freebayes_on_region
from pathlib import Path

TEMP_DIR = "buttermap_temp"
BAMS = ["test_data/test_output/sample1.MD.bam",  "test_data/test_output/sample2.MD.bam"]
REFERENCE = "test_data/test_input/brenthis_ino_chr14.reference.fasta"
REGION = "brenthis_ino.SP_BI_364.chromosome_14:0-100000"

Path(TEMP_DIR).mkdir(exist_ok=True)

def test_run_freebayes_on_region_no_zip_no_index():
    run_freebayes_on_region(bams = BAMS,
                            reference = REFERENCE,
                            region=REGION, 
                            output_dir=TEMP_DIR, bgzip=False, index=False)
    
    with open(f"{TEMP_DIR}/{REGION}.vcf") as f:
        line0 = f.readlines()[0]
    
    assert line0 == "##fileformat=VCFv4.2\n"

def test_run_freebayes_on_region_with_zip_with_index():
    run_freebayes_on_region(bams = BAMS,
                            reference = REFERENCE,
                            region=REGION, 
                            output_dir=TEMP_DIR, bgzip=True, index=True)
    
    assert Path(f"{TEMP_DIR}/{REGION}.vcf.gz").exists() is True


    
