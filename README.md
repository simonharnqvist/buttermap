# buttermap
Mapping and VCF calling pipeline for butterfly genomes from paired-end reads

## Pipeline information
### Pipeline steps
1. Read mapping with minimap2
2. BAM file generation with samtools
3. Duplicate marking with sambamba
4. Splitting VCF into regions
5. Variant calling per region with freebayes
6. Merging chromosome VCFs with bcftools

### Dependencies
|Software              | Version used in development |
|----------------------|-----------------------------|
|Python                | 3.11.5                      |
|ruffus                | 2.8.4                       |
|joblib                | 1.2.0                       |
|minimap2              | 2.26-r1175                  |
|samtools              | 1.17                        |
|sambamba              | 1.0.1                       |
|freebayes             | v1.3.6                      |
|bcftools              | 1.17                        |


### Input files

Buttermap requires the following input files, which **must** follow the naming rules below:

|File              | File name                   |
|------------------|-----------------------------|
|Reference FASTA   | [REFERENCE_SPECIES].reference.fasta   |
|Reads             | [SAMPLE_ID].[N].fastq.gz    |

## User guide

### Installation

1. Clone GitHub repo: </br> 
```git clone git+https://github.com/simonharnqvist/buttermap.git```

2. Install using conda:</br> 
```conda env create --file buttermap.yaml```

This should install `buttermap` itself, but if that doesn't work (or if you installed dependencies manually), try:
```pip install .``` in the top level directory.

Don't forget to `activate` the right conda environment after installation.

### Running `buttermap`
You should now be able to run `buttermap` from the command line. Before doing so, I would recommend creating a new directory, containing a subdirectory with **copies of** the input files. Some tasks in `buttermap` will change files in-place, and although the original input files *should* be untouched, there are no guarantees.


```
(buttermap) computer:dir buttermap
usage: buttermap [-h] [--reference-fasta REFERENCE_FASTA] [--reads-dir READS_DIR] [--temp-dir TEMP_DIR] [--threads THREADS] [--region-size REGION_SIZE]

Mapping and calling pipeline for butterfly genomes

options:
  -h, --help            show this help message and exit
  --reference-fasta REFERENCE_FASTA
                        Reference genome in FASTA format; filename must conform to [SPECIES].[...].fasta (default: None)
  --reads-dir READS_DIR
                        Directory containing FASTQ read files, each matching '[SAMPLE].[N].fastq.gz' format (default: None)
  --temp-dir TEMP_DIR   Directory for temporary files; automatically deleted after successful run (default: buttermap_temp)
  --threads THREADS     Number of CPU threads to use (default: 1)
  --region-size REGION_SIZE
                        Region size for parallelising variant calling (default: 100000)
```
### Output
`buttermap` will output a series of intermediary files in the same directory as the provided FASTQ files. It is a good idea to copy the FASTQ files to a new directory before running `buttermap`. 

The final output is a gzipped VCF `[REFERENCE_SPECIES].vcf.gz`, which is written to the current working directory unless piped elsewhere. Please ensure that this file does not already exist in your current working directory (in fact, make a new directory to run `buttermap` in, if possible), otherwise the pipeline will not run.
