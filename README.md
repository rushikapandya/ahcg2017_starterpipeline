# Applied Human Computational Biology 2017
## About me
I'm *Rushika Pandya* and I'm a second year MS Bioinformatics student. As a part of Applied Human Computational Genomics (AHCG) 2017, we are developing a variant calling pipeline for genomic data analysis.

## Mission
Our mission is to device liquid biopsy kit to for non-invasive early detection of cancer

Our project involves involves integrating wet lab and bioinformatics aspects of cancer detection. The pipeline that we have developed is for the bioinformatics analysis following the wet lab steps of exoDNA sequencing. We aim to optimize liquid biopsy variant calling pipelines making it more sensitive to rare variants and allow the detection and evaluation of somatic variants in early cancer stages. 

## ahcg_pipeline
Goal: Build variant calling pipeline for genomic data analysis to detect low frequency variants from liquid biopsy

The final version of the pipeline (v.1.0.8) is downloadable at [link](https://github.com/rushikapandya/ahcg2017_starterpipeline/blob/master/ahcg_pipeline.py)

### Features
* **Input:** Exome Sequence Data
* In-built auto retrieval of sample from SRA
* Variant filtering by quality and depth (QUAL>=30 and DP>=25)
* Calculates median, average, and max coverage of genes
* Detects copy number variants (CNVs) using ControlFreec
* **Output:** VCF file of variants & Summary result of CNVs

### Important Commands:

#### Build Directory Structure
```{sh}
mkdir -p data/reads data/reference data/adapters output
```

#### Run Pipeline
```{sh}
./ahcg_pipeline_v1.0.8.py -c config_file.txt
```
#### Help
```{sh}
./ahcg_pipeline.py -h
```

## Tool Requirements

1. [Python3 - version 3.4.1](https://www.python.org/download/releases/3.4.1/)
2. [Trimmomatic - version 0.36](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip)
3. [Bowtie2 - version 2.3.2](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.2)
4. [Picard tools - version 2.11.0](https://github.com/broadinstitute/picard/releases/download/2.11.0/picard.jar) (Note: Java Version 1.8 is required to run Picard)
5. [GATK - version 3.8](https://software.broadinstitute.org/gatk/download/)
6. [Samtools - version 1.6](https://downloads.sourceforge.net/project/samtools/samtools/1.6/samtools-1.6.tar.bz2?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fsamtools%2F&ts=1510018121&use_mirror=phoenixnap)
7. [Control-Freec version 11.0](https://github.com/BoevaLab/FREEC/archive/v11.0.tar.gz)
8. [R Language version 3.3.2](https://cran.cnr.berkeley.edu/) 

## Reference Genome
### Download Reference Genome
Reference genomes can be downloaded from [Illumina iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.html)
Source: NCBI
Build: GRCh38

### Create Index File
FASTA index file is created using Samtools
```{sh}
samtools faidx reference.fa
```
### Create Dictionary File
FASTA sequence dictionary file is created using Picard
```{sh}
java -jar picard.jar CreateSequenceDictionary R=reference O=dictionary
```

## Input data
### Test Data Used
* [SRR1654210 (tumor exome data)](https://www.ncbi.nlm.nih.gov/sra/?term=SRR1654210)
* [SRR1654222 (germline exome data)](https://www.ncbi.nlm.nih.gov/sra/SRR1654222)

For the input data, one can specify the path to the sample file or provide the SRA id in the config file

### Specify Path
Specify the path to sample SRA file using `inputfiles` option in the config file

### Auto Retrieve Sample from SRA
Specify the id of the sample SRA file (example: SRR948994) using the `sraid` option in the config file to delegate the download of the sample files to the pipeline

## Configuration file

The configuration file has two required sections `[tools]` and `[data]` in which the required software tools and input data are described. An additional optional section `[freec-control]` can be included to specify options for the `[control]` section of the configuration of the Control-FREEC tool. 

### `[data]` section

| Option       | Description                                                                  |
|--------------|------------------------------------------------------------------------------|
| `inputfiles` | List of paired end read files (comma sparated)                               |
| `geneset`    | Path to the bed file with genes of interest to calculate coverage statistics |
| `outputdir`  | Path to the output directory                                                 |
| `adapters`   | Path to adapters fasta file to perform sequence trimming with Trimmomatic    |
| `chrlenfile` | Path to file with chromosome lengths for Control-FREEC                       |
| `chrfiles`   | Path to the directory with chromosomes fasta files for Control-FREEC         |
| `dbsnp`      | Path to dbSNP vcf file for GATK                                              |
| `index`      | Path to the prefix of the reference Bowtie2 index                            |
| `reference`  | Path to the reference genome fasta file                                      |

### `[tools]` section

| Option        | Description                                                                                                |
|---------------|------------------------------------------------------------------------------------------------------------|
| `bowtie2`     | Path to Bowtie2 executable                                                                                 |
| `freec`       | Path to Control-FREEC executable                                                                           |
| `gatk`        | Path to GATK jar file                                                                                      |
| `makegraph`   | Path to Control-FREEC `makeGraph.R` script (Usually in the folder `scripts` at the Control-FREEC root dir) |
| `picard`      | Path to Picard jar file                                                                                    |
| `samtools`    | Path to Samtools executable                                                                                |
| `trimmomatic` | Path to Trimmomatic jar file                                                                               |

### `[freec-control]` section

Optional section for Control-FREEC's config file `control` section parameters

| Option          | Description                                                                                                                                                           |
| --------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| mateFile        | Path to file to act as control of the current sample. See control-FREEC manual for details>                                                                           |
| inputFormat     | Format of mateFile (SAM, BAM, pileup and others. See control-FREEC manual for details)>                                                                               |
| mateOrientation | Orientation of reads in mateFile. 0 - single ends), RF - Illumina mate-pairs, FR - Illumina paired-ends), FF - SOLiD mate-pairs. See control-FREEC manual for details |

### Example of config file

```
[data]
sraid           = SRR1654210
geneset         = /data2/AHCG2017FALL/guardant360/guardant360.refGene_hg38.genes.bed
outputdir       = /data2/AHCG2017FALL/output5

adapters        = /data2/AHCG2017FALL/bin/Trimmomatic-0.36/adapters/NexteraPE-PE.fa
chrlenfile      = /data2/AHCG2017FALL/reference_genome/chromosomeSizes.txt
chrfiles        = /data2/AHCG2017FALL/reference_genome/chroms/
dbsnp           = /data2/AHCG2017FALL/reference_genome/GATKResourceBundle/dbsnp_146.hg38.vcf.gz
index           = /data2/AHCG2017FALL/reference_genome/Bowtie2Index/genome
reference       = /data2/AHCG2017FALL/reference_genome/genome.fa

[tools]
assesssig       = /data2/AHCG2017FALL/bin/FREEC/scripts/assess_significance.R
bowtie2         = /data2/AHCG2017FALL/bin/bowtie2-2.2.9/bowtie2
fastq-dump      = /data2/AHCG2017FALL/bin/sratoolkit/bin/fastq-dump
freec           = /data2/AHCG2017FALL/bin/FREEC/src/freec
gatk            = /data2/AHCG2017FALL/bin/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar
java            = /data2/AHCG2017FALL/bin/java-1.8/bin/java
makegraph       = /data2/AHCG2017FALL/bin/FREEC/scripts/makeGraph.R
picard          = /data2/AHCG2017FALL/bin/picard/picard.jar
samtools        = /data2/AHCG2017FALL/bin/samtools-1.5/samtools
trimmomatic     = /data2/AHCG2017FALL/bin/Trimmomatic-0.36/trimmomatic-0.36.jar

[freec-control]
mateFile        = /data2/AHCG2017FALL/output4/SRR2530741_1_trimmed_final.bam
inputFormat     = BAM
mateOrientation = FR
```

## Filter Variants by Sequence Coverage
Post variant calling using the pipeline we filtered the variants based on depth of coverage

```{sh}
cat Variants.vcf | java -jar path/to/SnpSift.jar filter "((QUAL >= 30) && (DP >= 25))" > Variants-Filtered.vcf
```

## Additional Information

### Virtual Box Set-up
Click [here](https://www.perkin.org.uk/posts/create-virtualbox-vm-from-the-command-line.html) for detailed steps to install and set-up your virtual box

### Get Exon Co-ordinates
Click [here](https://github.com/rushikapandya/ahcg2017_starterpipeline/blob/master/get_exon_co-ordinates.pdf) for detailed instructions on downloading exon co-ordinates from UCSC genome browser


## Acknowledgments
Dr. Fredrik Vannberg and Cai Huang



