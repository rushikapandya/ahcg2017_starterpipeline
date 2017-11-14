# Applied Human Computational Biology 2017
## About me
I'm *Rushika Pandya* and I'm a second year MS Bioinformatics student. As a part of Applied Human Computational Genomics (AHCG) 2017, we are developing a variant calling pipeline for genomic data analysis.

## Mission
Our mission is to device liquid biopsy kit to for non-invasive early detection of cancer

Our project involves involves integrating wet lab and bioinformatics aspects of cancer detection. The pipeline that we have developed is for the bioinformatics analysis following the wet lab steps of exoDNA sequencing. We aim to optimize liquid biopsy variant calling pipelines making it more sensitive to rare variants. 

## ahcg_pipeline
Goal: Build variant calling pipeline for genomic data analysis to detect low frequency variants from liquid biopsy

## Filter Variants by Sequence Coverage
Post variant calling using the pipeline we filtered the variants based on depth of coverage

```{sh}
cat Variants.vcf | java -jar path/to/SnpSift.jar filter "((QUAL >= 30) && (DP >= 25))" > Variants-Filtered.vcf
```


## Requirements

1. [Python3 - version 3.4.1](https://www.python.org/download/releases/3.4.1/)
2. [Trimmomatic - version 0.36](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip)
3. [Bowtie2 - version 2.2.9](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/)
4. [Picard tools - version 2.6.0](https://github.com/broadinstitute/picard/releases/download/2.6.0/picard.jar)
5. [GATK - version 3.4](https://software.broadinstitute.org/gatk/download/)

## Reference genome

Reference genomes can be downloaded from [Illumina iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.html)

## Test data

Use the following protocol to download and prepare test dataset from NIST sample NA12878

```{sh}
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
head -100000 NIST7035_TAAGGCGA_L001_R1_001.fastq > test_r1.fastq
head -100000 NIST7035_TAAGGCGA_L001_R2_001.fastq > test_r2.fastq
```

## Help

To access help use the following command:

```{sh}
python3 ahcg_pipeline.py -h
```
