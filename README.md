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

Full list of dependencies are provided in `environment.yaml`. If you need to install dependencies manually (e.g. via homebrew for testing on Mac), please see this table:

| Software/package | Version |
| ---------------- | --------------------------- |
| Python           | 3.12.0                      |
| ruffus           | 2.8.4                       |
| joblib           | 1.3.2                       |
| minimap2         | 2.26                        |
| samtools         | 1.18                        |
| sambamba         | 1.0                       |
| freebayes        | 1.3.7                       |
| bcftools         | 1.18                        |

Additionally, pytest is required to run the unit tests, which is an optional installation step.


### Input files

Buttermap requires the following input files, which **must** follow the naming rules below:

| File            | File name                           |
| --------------- | ----------------------------------- |
| Reference FASTA | [REFERENCE_SPECIES].reference.fasta |
| Reads           | [SAMPLE_ID].[N].fastq.gz            |

## User guide

### Installation

#### 1. Clone GitHub repo: </br> 
``` shell
git clone https://github.com/simonharnqvist/buttermap.git
```

#### 2. Install using conda:</br> 
``` shell
cd buttermap # important; see below
conda env create --file environment.yaml # or use mamba
conda activate buttermap
```

**Important**: You must run this command from within the repo directory as it relies on a relative path to the `buttermap` wheel.

This should install `buttermap` itself, but if that doesn't work (or if you installed dependencies manually), try
```pip install .``` in the top level directory.

#### 3. (Optional): Run tests </br>
``` shell
python -m pytest tests
```
This is a good idea to ensure that the installation has worked and that the code wasn't broken to begin with. Please raise an issue on GitHub with the error message if the tests failed. 

Note that the tests are minimal and essentially just check that the pipeline doesn't throw any errors on a (tiny) test dataset, and that a correctly formatted VCF is produced. You should not assume that this means that the pipeline is bug-free.


### Running `buttermap`
You should now be able to run `buttermap` from the command line. Before doing so, I would recommend creating a new directory, containing a subdirectory with **copies of** the input files. Some tasks in `buttermap` will change files in-place, and although the original input files *should* be untouched, there are no guarantees.

```shell
cp my_original_input buttermap_copies
```

Parameters are available by running `buttermap`:
``` shell
(buttermap) computer:dir buttermap --help
```

```
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

The final output is a gzipped VCF `[REFERENCE_SPECIES].vcf.gz`, which is written to the current working directory unless piped elsewhere. Please ensure that this file does not already exist in your current working directory  otherwise the pipeline will not run.
