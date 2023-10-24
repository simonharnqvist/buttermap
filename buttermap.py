import subprocess
from ruffus import transform, suffix
from fasta_generate_regions import generate_regions

# get these from a YAML config file?
starting_files = [("example_data/brenthis_ino_chr14.reference.fasta",
                  "example_data/SO_BI_375.subsamp.1.fastq.gz",
                   "example_data/SO_BI_375.subsamp.2.fastq.gz")]

@transform(input=starting_files,
           filter=suffix(".reference.fasta"),
           output=".paf"
           )
def map_reads(infiles, outfile, threads=1):
    """Map reads to reference with minimap2

    Args:
        infiles (list): List of input files [(reference, reads1, reads2), ...]
        outfile (path): PAF file path
        threads (int, optional): Number of threads to run minimap with. Defaults to 1.
    """

    sample_id = Path(infiles[1]).stem.split(".")[0]
    assert sample_id == Path(infiles[2]).stem.split(".")[0]

    minimap_cmd = ["minimap2",
                   "-ax", "sr", "-t", "1", 
                   "-R", f'@RG\\tID:{sample_id}\\tSM:{sample_id}']
    minimap_cmd.extend(list(infiles))

    subprocess.run(minimap_cmd, stdout=open(outfile, "w"), shell=False, check=True)


@split(fasta_index_file, ".bed")
def generate_fasta_regions(fasta_index_file, region_beds, size, bed_dir):
    generate_regions(fasta_index_file=fasta_index_file,
                     size=size,
                     bed_dir=bed_dir)


@transform(map_reads,
           suffix(".paf"),
           ".bam")
def generate_bam(input_paf, output_bam, threads, temp_dir):
    """View & sort PAF file to generate BAM file

    Args:
        input_paf (path): Input PAF path
        output_bam (path): Output BAM path
        threads (int): Number of threads to use
        temp_bam (path): Temporary BAM file to write to during processing
    """
    view = subprocess.run(["samtools", "view", "-q", "1", "-b", input_paf], 
                          stdout=subprocess.PIPE, 
                          shell=True,
                          check=True)
    
    output_buff = open(output_bam, "w")
    subprocess.run(["samtools", "sort", 
                    f"-@{threads}", 
                    "-m2G", 
                    "-T", f"{temp_dir}/temp.bam"], 
                          stdin=view.stdout, 
                          stdout=output_buff, 
                          shell=True, check=True)
    
    subprocess.run(["samtools", "index", f"{output_bam}"], shell=True, check=True)


@transform(generate_bam,
           suffix(".bam"),
           ".MD.bam")
def mark_duplicates(input_bam, output_bam, threads, temp_dir, overflow_list_size=600_000):
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
        shell=True,
        check=True
    )


@split(mark_duplicates, "*.vcf", region_bed)
def call_variants_per_chromosome(input_bams_list, output_region_vcf, input_fasta, region):

    # How to parallelise this? For loop with joblib?
    
    variant_calling = subprocess.run([
        "freebayes", 
        "--fasta-reference", input_fasta,
        "--region", region,
        "--limit-coverage", "250", 
        "--use-best-n-alleles", "8", 
        "--no-population-priors", 
        "--use-mapping-quality", 
        "--ploidy" "2", 
        "--haplotype-length", "-1",
        "--bam_list", input_bams_list
        ],
        shell=True,
        check=True,
        stdin=
        stdout=open(output_vcf, "w"))


@merge(call_variants_per_chromosome, "*.vcf", "calls.vcf")
def merge_chromosome_vcfs(input_vcfs, output_vcf):
    subprocess.run(["bcftools", "concat", "--allow-overlaps", "--remove-duplicates"].append(input_vcfs),
                   stdout=output_vcf,
                   check=True,
                   shell=True)