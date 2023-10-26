#!/usr/bin/env python

from ruffus import transform, suffix, pipeline_run, merge, collate, regex, formatter
from pathlib import Path
import subprocess
from joblib import Parallel, delayed
import glob
import shutil
import argparse
import sys
from buttermap.fasta_generate_regions import generate_regions
from buttermap.utils import read_bed

#############
### SETUP ###
#############


parser = argparse.ArgumentParser(description="Mapping and calling pipeline for butterfly genomes",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--reference-fasta", 
                    help="Reference genome in FASTA format; filename must conform to [SPECIES].[...].fasta", type=str)
parser.add_argument("--reads-dir", 
                    help="Directory containing FASTQ read files, each matching '[SAMPLE].[N].fastq.gz' format")
parser.add_argument("--temp-dir", 
                    help="Directory for temporary files; automatically deleted after successful run", 
                    default="buttermap_temp")
parser.add_argument("--threads", 
                    help="Number of CPU threads to use",
                    default=1)
parser.add_argument("--region-size",
                    help="Region size for parallelising variant calling", default=100_000, required=False)
args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

REFERENCE = args.reference_fasta
READS_DIR = args.reads_dir
TEMP_DIR = args.temp_dir
THREADS = args.threads
REF_SPECIES = Path(REFERENCE).name.split(".")[0]
REGION_SIZE = args.region_size
READS = glob.glob(f"{READS_DIR}/*.fastq.gz")

Path(TEMP_DIR).mkdir(exist_ok=True)



#####################################################
### SUBPIPELINE 1: Index FASTA, generate regions file
#####################################################

@transform(REFERENCE, suffix(".fasta"), ".fasta.fai")
def index_fasta(infile, outfile):
    """Index FASTA file"""
    subprocess.run(["samtools", "faidx", infile, "--fai-idx", outfile], check=True)


@transform(index_fasta, suffix("fasta.fai"), "regions.bed", REGION_SIZE)
def generate_fasta_regions(infile, outfile, size):
    """Generate BED file of regions for parallelisation of downstream steps"""
    generate_regions(fasta_index_file=infile,
                     size=size, output_bed=outfile)


########################################
### SUBPIPELINE 2: Map reads per sample
########################################

@collate(READS, 
         formatter("([^/]+)[12].fastq.gz$"), 
         "{path[0]}/{1[0]}paf",
         REFERENCE, 1)
def map_reads(infiles, outfile, reference_fasta, threads):
    """Map reads to reference with minimap2

    Args:
        infiles (list): List of input files [reference, (sample1_reads1, sample1_reads2)...]
        outfile (path): PAF file path
        threads (int, optional): Number of threads to run minimap with. Defaults to 1.
    """

    sample_id = Path(infiles[0]).stem.split(".")[0]
    assert sample_id == Path(infiles[1]).stem.split(".")[0]

    minimap_cmd = ["minimap2",
                   "-ax", "sr", "-t", str(threads), 
                   "-R", f'@RG\\tID:{sample_id}\\tSM:{sample_id}']
    minimap_cmd.append(reference_fasta)
    minimap_cmd.extend(list(infiles))

    outbuff = open(outfile, "w")
    subprocess.run(minimap_cmd, check=True, stdout=outbuff)



@transform(map_reads, suffix(".paf"), ".bam", 1, TEMP_DIR)
def generate_bam(infile, outfile, threads, temp_dir):
    """View & sort PAF file to generate BAM file

    Args:
        input_paf (path): Input PAF path
        output_bam (path): Output BAM path
        threads (int): Number of threads to use
        temp_bam (path): Temporary BAM file to write to during processing
    """
    view = subprocess.Popen(["samtools", "view", "-q", "1", "-b", infile], 
                          stdout=subprocess.PIPE)
    
    subprocess.Popen(["samtools", "sort", 
                    f"-@{threads}", 
                    "-m2G", 
                    "-T", f"{temp_dir}/temp.bam"], 
                    stdin=view.stdout, 
                    stdout=open(outfile, "w"))
    
    #subprocess.run(["samtools", "index", outfile], check=True) # this fails - see if later steps require it


@transform(generate_bam, suffix(".bam"), ".MD.bam", 1, "temp", 600_000)
def mark_duplicates(input_bam, output_bam, threads, temp_dir, overflow_list_size):
    """Mark duplicates with Sambamba

    Args:
        input_bam (path): Input BAM file to mark duplicates on
        output_bam (path): Output BAM with marked duplicates
        threads (int): Number of threads to use
        temp_dir (path): Directory to write temporary output to
        overflow_list_size (int, optional): See Sambamba parameters. Defaults to 600_000.
    """
    subprocess.run(
        [
            "sambamba", "markdup", 
            "--overflow-list-size", f"{overflow_list_size}",
            "-t", f"{threads}",
            "--tmpdir", temp_dir,
            input_bam, output_bam
        ],
        check=True
    )


################################
### SUBPIPELINE 3: Call variants
################################

def run_freebayes_on_region(bams, reference, region, output_dir):
    """Run Freebayes on a single region. 
    Uses a non-Ruffusonic workaround (since no suitable decorator) to write VCFs to a temp dir outside a task.


    Args:
        bams (list): List of paths to BAM files
        reference (path): Path to reference FASTA
        region (str): Region in format {chrom}:{start}-{stop}
        output_dir (path): Directory to put region VCFs in
    """
    cmd = ["freebayes", "-f", reference, 
                    "--limit-coverage", "250",
                    "--use-best-n-alleles", "8",
                    "--no-population-priors", "--use-mapping-quality",
                    "--ploidy", "2", "--haplotype-length", "-1",
                    "--region", region]

    for bam in bams:
        cmd.extend(["--bam", bam])

    output_vcf = f"{output_dir}/{region}.vcf"

    subprocess.Popen(cmd, stdout=open(output_vcf, "w"))

    subprocess.run(["bgzip", "--force", output_vcf], check=True)
    subprocess.run(["bcftools", "index", f"{output_vcf}.gz"], check=True)


@merge([mark_duplicates, generate_fasta_regions], f"{REF_SPECIES}.vcf.gz", 
       REFERENCE, THREADS, TEMP_DIR)
def call_variants(infiles, outfile, reference, threads, tempdir):
    """Call variants with Freebayes.

    Args:
        infiles (list): Outputs of mark_duplicates and generate_fasta_region tasks
        outfile (_type_): Output VCF path}
        reference (_type_): Reference FASTA
        threads (_type_): Number of threads for parallel processing
        tempdir (_type_): Temporary directory to write VCFs to
    """

    Path(tempdir).mkdir(exist_ok=True)

    bams = infiles[:-1]
    regions = read_bed(infiles[-1])

    Parallel(n_jobs=threads)(delayed(run_freebayes_on_region)(bams, reference, region, tempdir) for region in regions[0:5])

    infiles = glob.glob(f"{tempdir}/*.vcf.gz")
    bcftools_concat_cmd = ["bcftools", "concat", "--allow-overlaps", "--remove-duplicates", "--output", outfile]
    bcftools_concat_cmd.extend(infiles)
    subprocess.run(bcftools_concat_cmd, check=True)

def main():

    Path(TEMP_DIR).mkdir(exist_ok=True)
    pipeline_run()
    shutil.rmtree(TEMP_DIR)

if __name__ == "__main__":
    main()