cadd2vcf
========

`cadd2vcf` is a tool to assist in using [CADD](https://github.com/kircherlab/CADD-scripts) scores to modify and annotate a VCF to VCF [specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf). It will also deal with situations where there are differences in contig chr-prefix between the VCF and CADD scores. More information on CADD scores can be found [here](https://cadd.gs.washington.edu/). 

## Installation
Installation is available via `conda`.
```
conda install -c bioconda cadd2vcf
```

## Usage
```
Usage: ./cadd2vcf [options] <vcf(can be gzipped)> <cadd.tsv (can be gzippped)>
cadd2vcf: Annotate vcf with CADD Report
-b --build Genome build used to generate CADD file (GRCh37 or GRCh38) default: GRCh38
-h  --help This help information.
```