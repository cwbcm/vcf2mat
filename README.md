# VCF to EA matrix
Repo compute EA matrix from VCF file. Code adapted from mu_calculator

## Installation
You can either use BigPipeline environment or follow these steps:
1. git clone https://github.com/cwbcm/vcf2mat.git
2. conda env create -f ./vcf2mat/environment.yml
3. conda activate vcf2mat


## Usage
Required arguments:
| Argument                | Descripion |
| ---------------------- |--------------------- |
| --VCF                | Path to annotated VCF file |
| --savepath           | Path for output file |
| --maxaf  | Sets maximum allele frequency threshold |
| --minaf  | Sets minimum allele frequency threshold |
| --cores              | Number of cpus to use |
| --output | Specify what EA matrices to output, use , to separate the outputs. E.g. sumEA,pEA,silent,maxEA |
| --format | Specify output file format (options: tsv, parquet). Parquet is recommended |

Optional arguments:
| Argument                 | Descripion |
| ---------------------- |--------------------- |
| --Ann      | Variant annotation pipeline used (options: ANNOVAR, VEP / default: VEP). Only VEP is currently supported |
| --ref      | Genome reference (options: hg19, hg38 / default: hg38) |
| --prefix   | Optional prefix to add before the EA matrices. If unspecified, sumEA_mat.parquet/tsv will be used |




## Command line example
```bash
#set your working directory to vcf2mat
cd ./vcf2mat
#run vcf2mat.py
python vcf2mat.py --VCF Path/to/vcf_file.vcf.gz --ref hg38 --savepath save/directory/ --cores 20 --maxaf 0.01 --format parquet
```

## Output
1. EA matrices (sumEA, pEA, maxEA or silent) for a given EA annotated vcf file. The genes are on the columns, while the samples are on the rows.
2. All samples in the vcf file will be included in the output. If only interested in a subset of samples, subset the vcf files first with bcftools.
3. Use custom reference files when only a subset of gene is of interest.


