import subprocess
from pathlib import Path

def run_freebayes_on_region(bams, reference, region, output_dir, bgzip=True, index=True):
    """Run Freebayes on a single region. 
    Uses a non-Ruffusonic workaround (since no suitable decorator) to write VCFs to a temp dir outside a task.


    Args:
        bams (list): List of paths to BAM files
        reference (path): Path to reference FASTA
        region (str): Region in format {chrom}:{start}-{stop}
        output_dir (path): Directory to put region VCFs in
    """

    Path(output_dir).mkdir(exist_ok=True)

    OUTPUT_VCF = f"{output_dir}/{region}.vcf"

    cmd = ["freebayes", "-f", reference, 
                    "--limit-coverage", "250",
                    "--use-best-n-alleles", "8",
                    "--no-population-priors", "--use-mapping-quality",
                    "--ploidy", "2", "--haplotype-length", "-1",
                    "--region", region,
                    "--vcf", OUTPUT_VCF]

    for bam in bams:
        cmd.extend(["--bam", bam])

    print(cmd)

    subprocess.run(cmd, check=True)

    if bgzip is True:
        subprocess.run(["bgzip", "--force", OUTPUT_VCF], check=True)
        OUTPUT_VCF = f"{OUTPUT_VCF}.gz"
    
    if index is True:
        subprocess.run(["bcftools", "index", OUTPUT_VCF], check=True)