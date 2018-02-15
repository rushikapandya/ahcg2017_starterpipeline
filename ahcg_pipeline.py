#!/usr/bin/env python3

import os
import sys
import glob
import logging
import argparse
import shutil
import subprocess
import statistics as st
import configparser as cp
import textwrap

# ==============================================================================
def calculateCoverageStats(bamFilePath, bedFilePath, outFilePath, samtoolsPath):
    """Calculate mean, median and average coverage of each region
    in bedFilePath by calling samtools depth, and save results to outFilePath."""
    covFile  = open(outFilePath, 'w')
    covFile.write("#chr\tstart\tstop\tname\tscore\tstrand\tmedian_cov\tavg_cov\tmax_cov\n")
    with open(bedFilePath) as fp:
        line = fp.readline()
        while line:
            tokens     = line.split("\t")
            region     = "{}:{}-{}".format(tokens[0],tokens[1],tokens[2])
            samcmd     = [samtoolsPath, "depth", "-r", region, bamFilePath]
            covLines   = subprocess.check_output(samcmd, universal_newlines=True).splitlines()
            perBaseCov = []
            for l in covLines:
                tokens = l.split("\t")
                perBaseCov.append(int(tokens[2]))

            covMedian  = 0
            covAverage = 0
            covMax     = 0

            if len(perBaseCov) > 0:
                covMedian  = st.median(perBaseCov)
                covAverage = st.mean(perBaseCov)
                covMax     = max(perBaseCov)

            covFile.write("{0}\t{1:.2f}\t{2:.2f}\t{3}\n".format(line.rstrip(),
                covMedian, covAverage, covMax))

            line = fp.readline()

    covFile.close()
# ==============================================================================
def getlist(option, sep=',', chars=None):
    """Return a list from a ConfigParser option. By default,
       split on a comma and strip whitespaces."""
    return [ chunk.strip(chars) for chunk in option.split(sep) ]
# ==============================================================================
def checkTool(key, name, confOptions):
    """Check if the executable file key defined in confOptions or present the
    system PATH, exists"""
    if key in confOptions:
        confOptions[key] = os.path.abspath(confOptions[key])
        if not os.path.exists(confOptions[key]):
            sys.stderr.write('ERROR: {} not found at {}\n'.format(name, confOptions[key]))
            sys.exit(1)
    else:
        toolPath = shutil.which(key)
        if toolPath != None:
                confOptions[key] = toolPath
        else:
            sys.stderr.write('ERROR: {} not found in system and not provided in the config file.\n'.format(name))
            sys.exit(1)
# ==============================================================================
def checkDataFile(key, name, confOptions):
    """Check if the data file defined by key in confOptions, exists"""
    if key in confOptions:
        confOptions[key] = os.path.abspath(confOptions[key])
        if key == 'index':
            indicies = glob.glob('{0}.*.bt2'.format(confOptions[key]))
            if len(indicies) == 0:
                sys.stderr.write('ERROR: {} not found at {}\n'.format(name, confOptions[key]))
                sys.exit(1)
        else:
            if not os.path.exists(confOptions[key]):
                sys.stderr.write('ERROR: {} not found at {}\n'.format(name, confOptions[key]))
                sys.exit(1)
    else:
        sys.stderr.write('ERROR: {} ({}) not found not provided in the config file\n'.format(name, key))
        sys.exit(1)
# ==============================================================================
def createControlFreecConfigFile(bamPath, confOptions, freecConfPath):
    '''Creates a config file for Control-FREEC based on files generated in this
    pipeline'''
    freecConf = cp.RawConfigParser()
    freecConf.optionxform = lambda option: option
    freecConf['general'] = {
        'chrLenFile' : confOptions['data']['chrlenfile'],
        'outputDir'  : '{}/freecOut'.format(os.path.abspath(confOptions['data']['outputdir'])),
        'chrFiles'   : confOptions['data']['chrfiles'],
        'samtools'   : confOptions['tools']['samtools'],
        'ploidy'     : 2,
        'window'     : 5000,
    }
    freecConf['sample'] = {
        'mateFile'        : bamPath,
        'inputFormat'     : 'BAM',
        'mateOrientation' : 'FR'
    }
    # if 'freec-sample' in confOptions:
    #     if 'mateOrientation' in confOptions['freec-sample']:
    #         freecConf['sample']['mateOrientation'] = confOptions['freec-sample']['mateOrientation']
    # else
    #     freecConf['sample']['mateOrientation'] = 0

    if 'freec-control' in confOptions:
        freecConf['control'] = confOptions['freec-control']

    with open(freecConfPath, 'w') as configfile:
        freecConf.write(configfile)
# ==============================================================================
def main(config_path):
    # Process config file
    config_path = os.path.abspath(config_path)
    config      = cp.ConfigParser(allow_no_value = True)
    config.read(config_path)
    confTools   = config['tools']
    confData    = config['data']
    toolsToCheck = {
        'assesssig'  : 'Control-FREEC\'s assess_significance.R script',
        'bowtie2'    : 'Bowtie2',
        'fastq-dump' : 'fastq-dump (SRA Toolkit)',
        'freec'      : 'Control-FREEC',
        'gatk'       : 'GATK',
        'java'       : 'Java',
        'makegraph'  : 'Control-FREEC\'s makeGraph.R script',
        'picard'     : 'Picard',
        'samtools'   : 'Samtools',
        'trimmomatic': 'Trimmomatic'
    }
    dataToCheck = {
        'adapters'   : 'Adapters file',
        'dbsnp'      : 'dbSNP vcf file',
        'geneset'    : 'Gene set bed file (genes of interest)',
        'index'      : 'Reference Bowtie index prefix',
        'reference'  : 'Reference genome file'
    }

    # ==========================================================================
    # Check for required tools
    for key in toolsToCheck:
        checkTool(key, toolsToCheck[key], confTools)

    # Check for required files (other than inputfiles, and outputdir options)
    for key in dataToCheck:
        checkDataFile(key, dataToCheck[key], confData)

    # Check for outputdir option
    if not 'outputdir' in confData:
        config['outputdir'] = './out'

    confData['outputdir'] = os.path.abspath(confData['outputdir'])

    # Create the output directory
    if not os.path.exists(confData['outputdir']):
        os.makedirs(confData['outputdir'])

    # Check for input files (inputfiles or sraid options)
    if 'inputfiles' in confData:
        files = getlist(confData['inputfiles'])
        # print("{}".format(getlist(confData['inputfiles'])))
        files = [os.path.abspath(fs) for fs in files]
        for f in files:
            if not os.path.exists(f):
                sys.stderr.write('ERROR: Fastq file not found at {}\n'.format(f))
                sys.exit(1)
    elif 'sraid' in confData:
        # Download SRA sample
        fqdumpCmd = [confTools['fastq-dump'], '--split-files', '-O', confData['outputdir'], confData['sraid']]
        print("INFO: Downloading SRA sample {}".format(confData['sraid']))
        fqDumpRun = subprocess.Popen(fqdumpCmd, shell=False)
        fqDumpRun.wait()

        if fqDumpRun.returncode != 0:
            sys.stderr.write('ERROR: SRA sample download failed; Exiting program\n')
            sys.exit(1)

        print("INFO: SRA sample {} downloaded to {}\n".format(confData['sraid'], confData['outputdir']))
    else:
        sys.stderr.write('ERROR: inputfiles or sraid options not provided in config file. Please specify one of them.\n')
        sys.exit(1)

    # ==========================================================================
    # Define fastq files depending on the source (sraid or inputfiles option)
    if 'inputfiles' in confData:
        inputFiles = getlist(confData['inputfiles'])
        read1      = inputFiles[0]
        read2      = inputFiles[1]
    else:
        read1 = '{}/{}_1.fastq'.format(confData['outputdir'], confData['sraid'])
        read2 = '{}/{}_2.fastq'.format(confData['outputdir'], confData['sraid'])

    # ==========================================================================
    # Trim fastq files
    tread1 = '{1}_trimmed.fq'.format(confData['outputdir'], os.path.splitext(read1)[0])
    tread2 = '{1}_trimmed.fq'.format(confData['outputdir'], os.path.splitext(read2)[0])
    sread1 = '{1}_unused.fq'.format(confData['outputdir'], os.path.splitext(read1)[0])
    sread2 = '{1}_unused.fq'.format(confData['outputdir'], os.path.splitext(read2)[0])

    tcmd   = [confTools['java'], '-jar', confTools['trimmomatic'], 'PE', '-phred33', read1, read2, tread1,
              sread1, tread2, sread2, 'ILLUMINACLIP:{0}:2:30:10'.format(confData['adapters']),
              'LEADING:0', 'TRAILING:0', 'SLIDINGWINDOW:4:15', 'MINLEN:36']

    print("INFO: Sequence reads trimming - Started")
    trun = subprocess.Popen(tcmd, shell=False)
    trun.wait()

    if trun.returncode != 0:
        sys.stderr.write('ERROR: Fastq trimming failed; Exiting program\n')
        sys.exit(1)

    print("INFO: Sequence reads trimming - Done!\n")

    # ==========================================================================
    # Align the reads using Bowtie2
    sam_path = '{1}.sam'.format(confData['outputdir'], os.path.splitext(tread1)[0])
    bcmd     = [ confTools['bowtie2'], '-x', confData['index'], '-S',
                 sam_path, '-p', '1' , '-1', tread1, '-2', tread2]

    print("INFO: Sequence read alignment - Started")
    brun = subprocess.Popen(bcmd, shell=False)
    brun.wait()

    if brun.returncode != 0:
        sys.stderr.write('ERROR: Bowtie2 failed; Exiting program\n')
        sys.exit(1)

    print("INFO: Sequence read alignment - Done!\n")
    # ==========================================================================
    # Add read group information
    add_path = '{0}/{1}_RG.bam'.format(confData['outputdir'],
        os.path.splitext(os.path.basename(sam_path))[0])
    acmd     = [confTools['java'], '-Xmx2g', '-jar', confTools['picard'],
                'AddOrReplaceReadGroups', 'I='+sam_path , 'O='+add_path,
                'SORT_ORDER=coordinate', 'RGID=Test', 'RGLB=ExomeSeq',
                'RGPL=Illumina', 'RGPU=HiSeq2500', 'RGSM=Test',
                'RGCN=AtlantaGenomeCenter', 'RGDS=ExomeSeq',
                'RGDT=2016-08-24', 'RGPI=null',
                'RGPG=Test', 'RGPM=Test', 'CREATE_INDEX=true']

    print("INFO: Read group information addition - Started")
    arun = subprocess.Popen(acmd, shell=False)
    arun.wait()

    if arun.returncode != 0:
        print('ERROR: Picard add read groups failed; Exiting program\n')
        sys.exit(1)

    print("INFO: Read group information addition - Done!\n")
    # ==========================================================================
    # Mark PCR duplicates
    dup_path = '{0}/{1}_MD.bam'.format(confData['outputdir'],
        os.path.splitext(os.path.basename(sam_path))[0])
    met_path = '{0}/{1}_MD.metrics'.format(confData['outputdir'],
        os.path.splitext(os.path.basename(sam_path))[0])
    mdcmd    = [confTools['java'], '-Xmx2g', '-jar', confTools['picard'], 'MarkDuplicates',
                'I='+add_path, 'O='+dup_path, 'METRICS_FILE='+met_path,
                'REMOVE_DUPLICATES=false', 'ASSUME_SORTED=true',
                'CREATE_INDEX=true']

    print("INFO: Duplicates marking - Started")
    mdrun = subprocess.Popen(mdcmd, shell=False)
    mdrun.wait()

    if mdrun.returncode != 0:
        print('ERROR: Picard mark duplicate failed; Exiting program\n')
        sys.exit(1)

    print("INFO: Duplicates marking - Done!\n")
    # ==========================================================================
    # Fix mate information
    fix_path = '{0}/{1}_FM.bam'.format(confData['outputdir'],
        os.path.splitext(os.path.basename(sam_path))[0])
    fcmd     = [confTools['java'], '-Xmx2g', '-jar', confTools['picard'],
                'FixMateInformation', 'I='+dup_path, 'O='+fix_path,
                'ASSUME_SORTED=true', 'ADD_MATE_CIGAR=true',
                'CREATE_INDEX=true']

    print("INFO: Mate information fixing - Started")
    frun = subprocess.Popen(fcmd, shell=False)
    frun.wait()

    if frun.returncode != 0:
        print('ERROR: Picard fix mate information failed; Exiting program\n')
        sys.exit(1)

    print("INFO: Mate information fixing - Done!\n")
    # ==========================================================================
    # Run realigner target creator
    interval_path = '{0}/{1}.intervals'.format(confData['outputdir'],
        os.path.splitext(os.path.basename(sam_path))[0])
    trcmd         = [confTools['java'], '-jar', confTools['gatk'],
                     '-T', 'RealignerTargetCreator', '-o', interval_path,
                     '-nt', '1', '-I', fix_path, '-R', confData['reference'],
                     '-known', confData['dbsnp']]

    print("INFO: Realignment targets creation - Started")
    trrun = subprocess.Popen(trcmd, shell=False)
    trrun.wait()

    if trrun.returncode != 0:
        print('Realigner Target creator failed; Exiting program\n')
        sys.exit(1)

    print("INFO: Realignment targets creation - Done!\n")
    # =========================================================================
    # Run indel realigner
    ral_path = '{0}/{1}_IR.bam'.format(confData['outputdir'],
        os.path.splitext(os.path.basename(sam_path))[0])
    recmd    = [confTools['java'], '-jar', confTools['gatk'], '-T', 'IndelRealigner',
                '--targetIntervals', interval_path, '-o', ral_path,
                '-I', fix_path, '-R', confData['reference']]

    print("INFO: Indel realignment - Started\n")
    rerun = subprocess.Popen(recmd, shell=False)
    rerun.wait()

    if rerun.returncode != 0:
        print('ERROR: Indel realigner creator failed; Exiting program\n')
        sys.exit(1)

    print("INFO: Indel realignment - Done!\n")
    # ==========================================================================
    # Base quality score recalibration
    bqs_path = '{0}/{1}.table'.format(confData['outputdir'],
        os.path.splitext(os.path.basename(sam_path))[0])
    bqscmd = [confTools['java'], '-jar', confTools['gatk'],
              '-T', 'BaseRecalibrator', '-R', confData['reference'],
              '-I', ral_path, '-o', bqs_path, '-nct', '1',
              '-cov', 'ReadGroupCovariate', '-knownSites', confData['dbsnp']]

    print("INFO: Base quality score recalibration - Started")
    bqsrun = subprocess.Popen(bqscmd, shell=False)
    bqsrun.wait()

    if bqsrun.returncode != 0:
        print('ERROR: Base quality score recalibrator failed; Exiting program\n')
        sys.exit(1)

    print("INFO: Base quality score recalibration - Done!\n")
    # ==========================================================================
    # Print Reads (generate final BAM)
    fbam_path = '{0}/{1}_final.bam'.format(confData['outputdir'],
        os.path.splitext(os.path.basename(sam_path))[0])
    prcmd     = [confTools['java'], '-jar', confTools['gatk'], '-T', 'PrintReads',
                 '-R', confData['reference'], '-I', ral_path,
                 '-o', fbam_path, '-BQSR', bqs_path, '-nct', '1']

    print("INFO: Final BAM generation strated")
    prrun = subprocess.Popen(prcmd, shell=False)
    prrun.wait()

    if prrun.returncode != 0:
        print('ERROR: Print reads failed; Exiting program\n')
        sys.exit(1)

    print("INFO: Final BAM generation - Done!\n")
    # ==========================================================================
    # Coverage stats calculation from the final BAM file
    covfile_path = '{0}/{1}_{2}_coverage.tsv'.format(
        confData['outputdir'],
        os.path.splitext(os.path.basename(confData['geneset']))[0],
        os.path.splitext(os.path.basename(fbam_path))[0])

    print("INFO: Coverage stats calculation - Started")
    calculateCoverageStats(fbam_path, confData['geneset'], covfile_path,
                           confTools['samtools'])
    print("INFO: Coverage stats calculation - Done!\n")

    # ==========================================================================
    # Haplotype caller
    vcf_path = '{0}/variants.vcf'.format(confData['outputdir'],
        os.path.splitext(os.path.basename(sam_path))[0])
    hcmd     = [confTools['java'], '-jar', confTools['gatk'],
                '-T', 'HaplotypeCaller', '-R', confData['reference'],
                '-I', fbam_path, '--dbsnp', confData['dbsnp'],
                '-o', vcf_path, '-nct', '1', '-gt_mode', 'DISCOVERY']

    print("INFO: Variants calling (Haplotype caller) - Started")
    hrun = subprocess.Popen(hcmd, shell=False)
    hrun.wait()

    if hrun.returncode != 0:
        print('ERROR: Haplotype caller failed; Exiting program\n')
        sys.exit(1)

    print("INFO: Variants calling (Haplotype caller) - Done!\n")
    # ==========================================================================
    # Variants filtering step
    vcf_filtered_path = '{0}/variants.filtered.vcf'.format(
        confData['outputdir'], os.path.splitext(os.path.basename(sam_path))[0])
    fcmd = [confTools['java'], '-jar', confTools['gatk'],
            '-T', 'SelectVariants', '-R', confData['reference'],
            '-V', vcf_path, '-select', 'DP>=25 && QUAL>=30',
            '-o', vcf_filtered_path]

    print("INFO: Variants filtering - Started")
    frun = subprocess.Popen(fcmd, shell=False)
    frun.wait()

    if frun.returncode != 0:
        print('ERROR: Variants filtering failed; Exiting program\n')
        sys.exit(1)

    print("INFO: Variants filtering - Done!\n")
    # ==========================================================================
    # CNVs calling
    freecConfPath = os.path.abspath('{0}/freec_conf.txt'.format(confData['outputdir']))

    createControlFreecConfigFile(os.path.abspath(fbam_path), config, freecConfPath)

    freecConf = cp.ConfigParser(allow_no_value = True)
    freecConf.read(freecConfPath)

    if not os.path.exists(freecConf['general']['outputDir']):
        os.mkdir(freecConf['general']['outputDir'])

    freecCmd = [confTools['freec'], '-conf', freecConfPath]

    print("INFO: CNVs calling - Started")
    freecRun = subprocess.Popen(freecCmd, shell=False)
    freecRun.wait()

    if freecRun.returncode != 0:
        print('ERROR: CNVs calling failed; Exiting program\n')
        sys.exit(1)

    print("INFO: CNVs calling - Done!\n")
    # ====================================================================
    # CNVs plots
    # cat /data2/AHCG2017FALL/bin/FREEC/scripts/makeGraph.R |
    # R --slave --args 2 outtest/freecOut/reads_1_trimmed_final.bam_ratio.txt
    ratioFilePath = os.path.abspath('{0}/{1}_ratio.txt'.format(freecConf['general']['outputDir'],
        os.path.basename(fbam_path)))
    catCmd  = ['cat', confTools['makegraph']]
    plotCmd = ['R', '--slave', '--args', freecConf['general']['ploidy'],
               ratioFilePath]

    print("INFO: CNVs plots generation - Started")
    catRun  = subprocess.Popen(catCmd, stdout=subprocess.PIPE, shell=False)
    plotRun = subprocess.Popen(plotCmd, stdin=catRun.stdout, shell=False)
    catRun.stdout.close()
    plotRun.wait()
    catRun.wait()

    if plotRun.returncode != 0:
        print('ERROR: CNVs plotting failed; Exiting program\n')
        sys.exit(1)

    print("INFO: CNVs plots generation - Done!\n")
    # ====================================================================
    # Significance assessing
    # cat /data2/AHCG2017FALL/bin/FREEC/scripts/assess_significance.R |
    # R --slave --args <*_CNVs> <*_ratio.txt>
    cnvsFilePath  = os.path.abspath('{0}/{1}_CNVs'.format(freecConf['general']['outputDir'],
        os.path.basename(fbam_path)))
    catCmd = ['cat', confTools['assesssig']]
    sigCmd = ['R', '--slave', '--args', cnvsFilePath, ratioFilePath]

    print("INFO: CNVs significance assessing - Started")
    catRun = subprocess.Popen(catCmd, stdout=subprocess.PIPE, shell=False)
    sigRun = subprocess.Popen(sigCmd, stdin=catRun.stdout, shell=False)
    catRun.stdout.close()
    sigRun.wait()
    catRun.wait()

    if sigRun.returncode != 0:
        print('ERROR: CNVs significance assessing failed; Exiting program\n')
        sys.exit(1)

    print("INFO: CNVs significance assessing - Done!\n")
    # ====================================================================

    print('INFO: Variant call pipeline completed')
    print('INFO: VCF file can be found at {0}'.format(vcf_path))
    return 0;

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        epilog = textwrap.dedent('''
config file format and options:
  [data]
  inputfiles  = <List of paired end read files (comma sparated) of the sample to process>
  sraid       = <SRA accession number of the sample to download and process>
  geneset     = <Path to the bed file with genes of interest to calculate coverage statistics>
  outputdir   = <Path to the output directory>

  adapters    = <Path to adapters fasta file to perform sequence trimming with Trimmomatic>
  chrlenfile  = <Path to file with chromosome lengths for Control-FREEC>
  chrfiles    = <Path to the directory with chromosomes fasta files for Control-FREEC>
  dbsnp       = <Path to dbSNP vcf file for GATK>
  index       = <Path to the prefix of the reference Bowtie2 index>
  reference   = <Path to the reference genome fasta file>

  [tools]
  assesssig   = <Path to Control-FREEC assess_significance.R>
                (Usually in the folder scripts at the Control-FREEC root dir)
  bowtie2     = <Path to Bowtie2 executable>
  freec       = <Path to Control-FREEC executable>
  fastq-dump  = <Path to fastq-dump (from SRA toolkit)>
  gatk        = <Path to GATK jar file>
  java        = <Path to Java executable>
  makegraph   = <Path to Control-FREEC makeGraph.R script>
                (Usually in the folder scripts at the Control-FREEC root dir)
  picard      = <Path to Picard jar file>
  samtools    = <Path to Samtools executable>
  trimmomatic = <Path to Trimmomatic jar file>

  [freec-control] # optional section for Control-FREEC app
  mateFile        = <Path to file to act as control of the current sample>
                    See control-FREEC manual for details.
  inputFormat     = <Format of mateFile> 
                    SAM, BAM, pileup and others. See control-FREEC manual 
                    for details.
  mateOrientation = <Orientation of reads in mateFile.> 
                    0 - single ends, RF - Illumina mate-pairs, 
                    FR - Illumina paired-ends, FF - SOLiD mate-pairs.
                    See control-FREEC manual for details.


config file details:
  - Sections [tools] and [data] are required
  - All the options in [tools] and [data] are required except for those that
    correspond to executable files (bowtie2, fastq-dump, freec, java, and
    samtools) which can be omitted if they are available in the system through
    the $PATH variable.
  - inputfiles and sraid options ([data] section) are mutually exclusive.
    If the two options are specified, inputfiles will have priority
'''))
    parser.add_argument('-c', '--config', dest = 'config_path', type = str,
        help = 'Path to config file')

    parser.add_argument('-v', '--version', action = 'version',
        version = '''%(prog)s 1.0.8
Copyright (c) 2017 Georgia Institute of Technology
Applied Human Computational Genomics - Fall 2017

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see http://www.gnu.org/licenses/.''')

    # Enrichment kit (eg Nextera Enrichment Capture kit)
    # Download samples given an SRA accession - Done!

    if len(sys.argv) == 1:
        parser.print_help()
        print('--------------------------------------------------------------------------------')
        print('{0} - ERROR: One or more required arguments are missing.\n'.format(parser.prog))
        sys.exit(1)

    args = parser.parse_args()

    sys.exit(main(args.config_path))
