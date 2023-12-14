#!/usr/bin/env python

from ruffus import transform, suffix, pipeline_run, merge, collate, regex, formatter, split, follows
from pathlib import Path
import subprocess
import glob
import shutil
import argparse
from math import ceil
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
parser.add_argument("--region-size",
                    help="Region size for parallelising variant calling", type=int, default=100_000, required=False)
parser.add_argument("--threads", help="Number of threads", default=1)
parser.add_argument("--output-dir",
                    help="Output directory", default=".")
args = parser.parse_args()

REFERENCE = args.reference_fasta
READS_DIR = Path(args.reads_dir)
TEMP_DIR = Path(args.temp_dir)
REF_SPECIES = Path(REFERENCE).name.split(".")[0]
REGION_SIZE = args.region_size
READS = glob.glob(f"{READS_DIR}/*.fastq.gz")
THREADS = args.threads
OUTPUT_DIR = args.output_dir



#####################################################
### SUBPIPELINE 1: Index FASTA, generate regions file
#####################################################

@transform(REFERENCE, suffix(".fasta"), ".fasta.fai")
def index_fasta(infile, outfile):
    """Index FASTA file"""
    subprocess.run(["samtools", "faidx", infile, "--fai-idx", outfile], check=True)


@split(index_fasta, f"{TEMP_DIR}/*.bed", REGION_SIZE)
def generate_fasta_regions(infile, outfiles, region_size):
    """Generate BED file per region for parallelisation of downstream steps. 
    Modified from https://github.com/freebayes/freebayes/blob/master/scripts/fasta_generate_regions.py"""


    with open(infile, "r") as fasta_index:
        for line in fasta_index:
            fields = line.strip().split("\t")
            chrom_name = fields[0]
            chrom_length = int(fields[1])

            region_start = 0
            while region_start < chrom_length:
                region_end = region_start + region_size
                if region_end > chrom_length:
                    region_end = chrom_length
                start = str(region_start)
                end = str(region_end)
                
                region = str(ceil(region_end / region_size))
                file_path = f"{TEMP_DIR}/{chrom_name}.{region}.bed"
                
                with open(file_path, "w") as f:
                    f.write("\t".join([chrom_name, start, end]))

                region_start = region_end


########################################
### SUBPIPELINE 2: Map reads per sample
########################################

@collate(READS, 
         formatter("([^/]+)[12].fastq.gz$"), 
         OUTPUT_DIR + "/{1[0]}sam",
         #"{path[0]}/{1[0]}sam",
         REFERENCE, THREADS)
def map_reads(infiles, outfile, reference_fasta, threads):
    """Map reads to reference with bwa-mem

    Args:
        infiles (list): List of input files [reference, (sample1_reads1, sample1_reads2)...]
        outfile (path): PAF file path
        threads (int, optional): Number of threads to run minimap with. Defaults to 1.
    """

    sample_id = Path(infiles[0]).stem.split(".")[0]
    assert sample_id == Path(infiles[1]).stem.split(".")[0]

    bwa_cmd = ["bwa", "mem", "-t", str(threads), "-R", f"@RG\\tID:{sample_id}\\tSM:{sample_id}"]
    bwa_cmd.append(reference_fasta)
    bwa_cmd.extend(list(infiles))

    print(bwa_cmd)

    outbuff = open(outfile, "w")
    subprocess.run(bwa_cmd, check=True, stdout=outbuff)


@transform(map_reads, suffix(".sam"), ".bam", THREADS, TEMP_DIR)
def generate_bam(infile, outfile, threads, temp_dir):
    """View & sort SAM file to generate BAM file

    Args:
        in (path): Input PAF path
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


@transform(generate_bam, suffix(".bam"), ".bam.bai")
def index_bam(infile, outfile):
    subprocess.run(["samtools", "index", infile], check=True)


@follows(index_bam)
@transform(generate_bam, suffix(".bam"), ".MD.bam", TEMP_DIR, 600_000, THREADS)
def mark_duplicates(input_bam, output_bam, temp_dir, overflow_list_size, threads):
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

#@follows(mark_duplicates)
@transform([generate_fasta_regions, mark_duplicates], suffix(".bed"), ".vcf")
def call_variants_on_region(infile, outfile):
    """Call variants on each region (i.e. from each BED)"""

    region = read_bed(infile)[0]
    bams = glob.glob(f"{OUTPUT_DIR}/*.MD.bam")

    cmd = ["freebayes", "-f", REFERENCE, 
                "--limit-coverage", "250",
                "--use-best-n-alleles", "8",
                "--no-population-priors", "--use-mapping-quality",
                "--ploidy", "2", "--haplotype-length", "-1",
                "--region", region,
                "--vcf", outfile]

    assert len(bams) >= 1, f"No BAMs found in {bams}"
    for bam in bams:
        cmd.extend(["--bam", bam])

    subprocess.run(cmd, check=True)


@transform(call_variants_on_region, suffix(".vcf"), ".vcf.gz")
def compress_vcf(infile, outfile):
    subprocess.Popen(["bgzip", "--force", infile])

@transform(compress_vcf, suffix(".vcf.gz"), ".vcf.gz.csi")
def index_vcf(infile, outfile):
    subprocess.Popen(["bcftools", "index", infile])

@follows(index_vcf)
@merge(compress_vcf, f"{OUTPUT_DIR}/{REF_SPECIES}.vcf.gz")
def merge_variant_calls(infiles, outfile):
    """Merge all region VCFs to single output"""

    bcftools_concat_cmd = ["bcftools", "concat", "--allow-overlaps", "--remove-duplicates", "--output", outfile]
    bcftools_concat_cmd.extend(infiles)
    subprocess.run(bcftools_concat_cmd, check=True, stdout=open(outfile, "w"))


##########
### MAIN 
##########

def main():
    Path(TEMP_DIR).mkdir()

    pipeline_run(multithread=THREADS)
    shutil.rmtree(TEMP_DIR)

if __name__ == "__main__":
    main()