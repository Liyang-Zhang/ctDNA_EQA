#!/bin/bash

# This script input tumor and cancer raw sequencing data (fastq format) and ouput pre-processed clean bam data for downstream snv calling
# It follows the guide from GATK best practise https://github.com/gatk-workflows/gatk4-data-processing/blob/master/processing-for-variant-discovery-gatk4.wdl

#SBATCH --mail-type=end
#SBATCH --mail-user=zly18810602991@163.com
#SBATCH --job-name=test
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH -J Preprocess
#SBATCH --ntasks-per-node=8
#SBATCH --output=%j.o
#SBATCH --error=%j.e

# Initialize a variable with an intuitive name to store the name of the input fastq file

R1=$1
R2=$2
sample=$3
project_dir=$4
pl=ILLUMINA


# Data from Pfizer
#baseSuper=`basename ${sample} .fastq`
baseSuper=`basename ${sample} .fastq.gz`
base=$(echo ${baseSuper} | cut -d'_' -f 1 | cut -d'-' -f 4)  # get the filename with basename function
echo ${base}

# Directory with reference genome
genome_19=/dssg/home/acct-medkwf/medkwf4/database/GATK/hg19/ucsc.hg19.fasta

# Set working dir
#working_dir=/dssg/home/acct-medkwf/medkwf4/results/MRD/CC_data/${base}
#working_dir=/dssg/home/acct-medkwf/medkwf4/results/MRD/huashanHospital/${base}
working_dir=/dssg/home/acct-medkwf/medkwf4/results/ctDNA_EQA/${base}
mkdir -p ${working_dir}

# Make dirs
mkdir -p ${working_dir}/FastQC/
mkdir -p ${working_dir}/bam_19/
mkdir -p ${working_dir}/interfiles/
mkdir -p ${working_dir}/gatk/
mkdir -p ${working_dir}/trim_galore

# Set up output file names
fastqc_out=${working_dir}/FastQC
bam=${working_dir}/bam_19
inter=${working_dir}/interfiles
gatk=${working_dir}/gatk
trim_galore_out=${working_dir}/trim_galore

# Set up file names
unmap_bamfile=${bam}/${base}.bam
intermediate_fq=${inter}/${base}.fq
aligned_bam=${bam}/${base}.align.bam
merge_bam=${bam}/${base}.merged.bam
tagfixed_bam=${bam}/${base}.tagfix.bam
dupmarked_bam=${bam}/${base}.dedup.bam
metrics_txt=${inter}/${base}.dedup.txt
recalibration_tab=${gatk}/${base}.bqsr.tab
bqsr_bam=${bam}/${base}.bqsr.bam
bqsr_sorted=${bam}/${base}.sort.bam

# Set up the software environment
module load miniconda3
source activate dna
#export IMAGE_NAME=/dssg/home/acct-medkwf/medkwf4/software/gatk.v4.2.5.0.sif
#module load picard
#module load bwa
#module load samtools

# specify the number of cores to use
cores=8

#echo "Processing file ${R1} ${R2}"
# Run FastQC
fastqc -t ${cores} \
${R1} \
-o ${fastqc_out}

fastqc -t ${cores} \
${R2} \
-o ${fastqc_out}

# Run trim_galore for paired end data #--fastqc \

trim_galore --quality 20 --phred33 \
--length 20 --stringency 1 \
-o ${trim_galore_out} \
--paired ${R1} ${R2} \
2> ${trim_galore_out}/${base}.trimgalore.log

# Convert fastq to unmapped bam files
# Make sure the input file names are correct
#picard FastqToSam F1=${trim_galore_out}/${base}_1_val_1.fq.gz F2=${trim_galore_out}/${base}_2_val_2.fq.gz O=${unmap_bamfile} SM=${sample} PL=${pl}
picard FastqToSam F1=${trim_galore_out}/${baseSuper}_val_1.fq.gz F2=${trim_galore_out}/${baseSuper/R1/R2}_val_2.fq.gz O=${unmap_bamfile} SM=${sample} PL=${pl}


# Convert bam to fastq and use bwa to align
# Offical script not work on server
picard SamToFastq I=${unmap_bamfile} FASTQ=${intermediate_fq}  INTERLEAVE=true NON_PF=true
bwa mem -K 100000000 -p -v 3 -t 16 -Y ${genome_19} ${intermediate_fq} | samtools view -1 - -o ${aligned_bam}
rm -fr ${intermediate_fq}

# Merge unmapped bam with mapped bam files
picard MergeBamAlignment \
VALIDATION_STRINGENCY=SILENT \
EXPECTED_ORIENTATIONS=FR \
ATTRIBUTES_TO_RETAIN=X0 \
ALIGNED_BAM=${aligned_bam} \
UNMAPPED_BAM=${unmap_bamfile} \
OUTPUT=${merge_bam} \
REFERENCE_SEQUENCE=${genome_19} \
PAIRED_RUN=true \
SORT_ORDER="unsorted" \
IS_BISULFITE_SEQUENCE=false \
ALIGNED_READS_ONLY=false \
CLIP_ADAPTERS=false \
MAX_RECORDS_IN_RAM=2000000 \
ADD_MATE_CIGAR=true \
MAX_INSERTIONS_OR_DELETIONS=-1 \
PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
PROGRAM_RECORD_ID="bwamem" \
PROGRAM_GROUP_VERSION="GetBwaVersion.version" \
PROGRAM_GROUP_COMMAND_LINE="bwa mem -K 100000000 -p -v 3 -t 16 -Y ${genome_19}" \
PROGRAM_GROUP_NAME="bwamem" \
UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
ALIGNER_PROPER_PAIR_FLAGS=true \
UNMAP_CONTAMINANT_READS=true

# Sort BAM file by coordinate order and fix tag values for NM and UQ
# This may be needed when MergeBamAlignment was run with SORT_ORDER other than 'coordinate' and thus could not fix these tags then
picard SortSam \
INPUT=${merge_bam} \
OUTPUT=/dev/stdout \
SORT_ORDER="coordinate" \
CREATE_INDEX=false \
CREATE_MD5_FILE=false \
| \
picard SetNmMdAndUqTags \
INPUT=/dev/stdin \
OUTPUT=${tagfixed_bam} \
CREATE_INDEX=true \
CREATE_MD5_FILE=false \
REFERENCE_SEQUENCE=${genome_19}

# Mark duplicate reads to avoid counting non-independent observations
# Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly.
# This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
# While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
# The maximum offset between two duplicate clusters in order to consider them optical duplicates. Default value 100
picard MarkDuplicates \
INPUT=${tagfixed_bam} \
OUTPUT=${dupmarked_bam} \
METRICS_FILE=${metrics_txt} \
VALIDATION_STRINGENCY=SILENT \
OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
ASSUME_SORT_ORDER="queryname" \
CREATE_MD5_FILE=false

# Generate Base Quality Score Recalibration (BQSR) model
# dbSNP_vcf and known_indels_sites_VCFs? where to get
# -L ~{sep=" -L " sequence_group_interval}
gatk --java-options "-Xmx128G" \
BaseRecalibrator \
-R ${genome_19} \
-I ${dupmarked_bam} \
--use-original-qualities \
-O ${recalibration_tab} \
--known-sites /dssg/home/acct-medkwf/medkwf4/database/GATK/hg19/1000G_omni2.5.hg19.sites.vcf.gz \
--known-sites /dssg/home/acct-medkwf/medkwf4/database/GATK/hg19/1000G_phase1.indels.hg19.sites.vcf.gz \
--known-sites /dssg/home/acct-medkwf/medkwf4/database/GATK/hg19/dbsnp_138.hg19.excluding_sites_after_129.vcf.gz \
--known-sites /dssg/home/acct-medkwf/medkwf4/database/GATK/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz \
--known-sites /dssg/home/acct-medkwf/medkwf4/database/GATK/hg19/dbsnp_138.hg19.vcf.gz \
--known-sites /dssg/home/acct-medkwf/medkwf4/database/GATK/hg19/hapmap_3.3.hg19.sites.vcf.gz \
--known-sites /dssg/home/acct-medkwf/medkwf4/database/GATK/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz

# Apply Base Quality Score Recalibration (BQSR) model
#-L ~{sep=" -L " sequence_group_interval} \
#singularity exec $IMAGE_NAME gatk --java-options "-Xmx128G" \
gatk --java-options "-Xmx128G" \
ApplyBQSR \
-R ${genome_19} \
-I ${dupmarked_bam} \
-O ${bqsr_bam} \
-bqsr ${recalibration_tab} \
--static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
--add-output-sam-program-record \
-OQ

# Index the bqsr.bam files
samtools index ${bqsr_bam}
samtools sort -o ${bqsr_sorted} -@ $cores ${bqsr_bam}
samtools index -b ${bqsr_sorted}

# Remove some intermediate files to save space
#rm -fr ${trim_galore_out}
#rm -fr ${unmap_bamfile}
#rm -fr ${aligned_bam}
#rm -fr ${merge_bam}
#rm -fr ${tagfixed_bam}

