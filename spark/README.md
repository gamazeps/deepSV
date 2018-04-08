# Read extraction process

This script is used for extracting the reads associated to structural variants.

It iterates through a given VCF file and extract the reads around the breakpoints of each variant it
contains.  

Regarding the data it assumes the following:

- vcf, aligned reads, and reference are using the same coordinates.
- The aligned reads for each sample are in a folder named after the sample, following the same
  structure as the one on the [1000 genome project ftp](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/).
- There is a *single* bam file in that folder (they are grabbed with `$SAMPLE/*bam`).
- The sample name is contained in the VCF INFO field (under SAMPLE).

## Configuration

The configuration file is written in json, the fields are the following:

- `supporting_reads_path`: folder in which we write the data.
- `reference_path`: path to the reference genome (in `fasta`)

## Fetching the data

In our use case we use the following sources for the data:

- Reference genome: GRch37 from 1000GP with `wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta`
- Variant calls from dbVar with `wget ftp://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/vcf/estd219.GRCh37.variant_call.vcf.gz`
- Bam files from 1000GP phase 3
