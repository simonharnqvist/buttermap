import subprocess
from ruffus import transform

@transform(input_fastas,
           suffix(".fasta"),
           ".paf")
def map_reads(input_fastas, output_paf, sample_id, threads):
    """Map reads to reference with minimap2

    Args:
        fastas (list): List [reference_fasta, read_fastas] of files
        output (path): Output path
        sample_id (str): Sample ID, e.g. "ES_BI_375"
        threads (int): Number of threads to use
    """

    minimap_cmd = ["minimap2", "-ax", "sr", 
                   "-t" f"{threads}" 
                   "-R" f"'@RG\tID:{sample_id}\tSM:{sample_id}'"]
    minimap_cmd.extend(input_fastas)
    minimap_cmd.append(f"> {output_paf}")

    subprocess.run(minimap_cmd, shell=True, check=True)


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

@transform(mark_duplicates,
           suffix(".MD.bam"),
           ".vcf")
def call_variants(input_files, output_vcf, region_size=100_000_000):
    ## TODO: make call_variants_per_region function and use this to paralellise
    """Call variants with Freebayes

    Args:
        input_files (list): List of input files [input_fasta, input_bams_list], where input_bams_list is a text file of BAMs
        output_vcf (path): Path to write VCF to
        region_size (int, optional): _description_. Defaults to 100_000_000.
    """

    input_fasta = input_files[0]
    input_bams_list = input_files[1]

    generate_regions = subprocess.run([
        "fasta_generate_regions.py", 
        "--fai", input_fasta,
        "--region_size", f"{region_size}"],
        stdout=subprocess.PIPE,
        shell=True,
        check=True)
    
    # Need to parallelise this - better to do this directly rather than using freebayes-parallel which isn't very robust
    variant_calling = subprocess.run([
        "freebayes", 
        "--fasta-reference", input_fasta,
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
