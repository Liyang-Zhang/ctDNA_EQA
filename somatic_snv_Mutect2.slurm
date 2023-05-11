#!/bin/bash

#SBATCH --mail-type=end
#SBATCH --mail-user=zly18810602991@163.com
#SBATCH --job-name=Mutect2
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH -J Mutect2
#SBATCH --ntasks-per-node=24
#SBATCH --output=%j.o
#SBATCH --error=%j.e

# Usage: sbatch somatic_snv_Mutect2.slurm /dssg/home/acct-medkwf/medkwf4/results/MRD/CC_data/SRR3862160_CC-H007C/bam_19/SRR3862160_CC-H007C.sort.bam /dssg/home/acct-medkwf/medkwf4/results/MRD/CC_data/SRR3862162_CC-H007N/bam_19/SRR3862162_CC-H007N.sort.bam

# sbatch somatic_snv_Mutect2.slurm /dssg/home/acct-medkwf/medkwf4/results/MRD/CC_data/SRR3862151_CC-H016C/bam_19/SRR3862151_CC-H016C.sort.bam /dssg/home/acct-medkwf/medkwf4/results/MRD/CC_data/SRR3862147_CC-H016N/bam_19/SRR3862147_CC-H016N.sort.bam

# Initialize a variable with an intuitive name to store the name of the input fastq file
tumor=$1
normal=$2
#project_dir=$3

# Get base name
baseSuper=`basename ${tumor} .sort.bam`
#base=$(echo ${baseSuper} | cut -d'_' -f 2)  # get the filename with basename function
base=${baseSuper}
echo "base is ${base}"
echo "Tumor sample is ${tumor}"
echo "Normal sample is ${normal}"

# Set working dir and annotation dir
#project_dir=/dssg/home/acct-medkwf/medkwf4/results/MRD/huashanHospital/${base}
#project_dir=/dssg/home/acct-medkwf/medkwf4/results/MRD/CC_data/${base}
project_dir=/dssg/home/acct-medkwf/medkwf4/results/ctDNA_EQA/${base}
mkdir -p ${project_dir}
working_dir=${project_dir}/Muetct2
mkdir -p ${working_dir}
annovar=${project_dir}/ANNOVAR
mkdir -p ${annovar}
interval_folder=${working_dir}/interval_files
mkdir -p ${interval_folder}
filter_folder=${working_dir}/filter
mkdir -p ${filter_folder}

# Set interval list file
#interval_list=${working_dir}/interval.list
#rm -fr ${interval_list}
#touch ${interval_list}


# Directory with reference genome
genome_19=/dssg/home/acct-medkwf/medkwf4/database/GATK/hg19/ucsc.hg19.fasta
#interval_list_bed=/dssg/home/acct-medkwf/medkwf4/database/GATK/hg19/ucsc.hg19.fasta.bed
#interval_list_bed=/dssg/home/acct-medkwf/medkwf4/reference/bed/UCSC_hg19_knowngene.bed
interval_list_bed=/dssg/home/acct-medkwf/medkwf4/reference/bed/602panel/NanOnco_Plus_Panel_v2.0_Covered_b37_cg.parY2X.sort.chr.bed

# Set up environment
module load miniconda3
source activate dna

# specify the number of cores to use
cores=24

# Set output files
normal_name=${working_dir}/${base}.normalname.txt
output_Mutect=${working_dir}/${base}.mutect2.vcf
output_Mutect_gz=${working_dir}/${base}.mutect2.vcf.gz
filter_Mutect=${working_dir}/${base}.filter.vcf
filter_Mutect_Man=${working_dir}/${base}.filter.mannually.vcf
annovar_aviinput_Mutect2=${annovar}/${base}.mutect2
annovar_annotated_Mutect2=${annovar}/${base}.mutect2.csv

# Use Mutect2
gatk --java-options "-Xmx128G" \
GetSampleName \
-R ${genome_19} \
-I ${normal} \
-O ${normal_name} \
-encode

# Directly run Mutect2 and fillter step
gatk --java-options "-Xmx128G" Mutect2 \
--native-pair-hmm-threads 24 \
-R ${genome_19} \
-I ${tumor} \
-I ${normal} \
-normal `cat ${normal_name}` \
-L ${interval_list_bed} \
-O ${output_Mutect_gz}



gatk --java-options "-Xmx128G" \
FilterMutectCalls \
-R ${genome_19} \
-V ${output_Mutect_gz} \
-O ${filter_Mutect} \
--min-median-base-quality 20 \
--min-median-mapping-quality 20


#########################################################################################
# split intervals 拆分参考基因组并不能节省时间

#test for spliting intervals to save running time
#gatk --java-options "-Xmx128G" \
#SplitIntervals \
#-R ${genome_19} \
#-L ${interval_list_bed} \
#--scatter-count 20 \
#-O ${interval_folder}
#
#
#for interval in `ls -l ${interval_folder} | grep .interval_list | awk '{print  $9}' | uniq`
#do
#  echo "interval file is ${interval_folder}/${interval}"
#  baseInterval=$(echo ${interval} | cut -d'_' -f 1)
#  echo "${output_Mutect_gz}.${baseInterval}" >> ${interval_list}
#  gatk --java-options "-Xmx128G" Mutect2 \
#  --native-pair-hmm-threads 24 \
#  -R ${genome_19} \
#  -I ${tumor} \
#  -I ${normal} \
#  -normal `cat ${normal_name}` \
#  -L ${interval_folder}/${interval} \
#  -O ${output_Mutect}.${baseInterval} &
#done
#wait
#echo "Finish parallelly running Mutect2"
#
#while IFS="" read -r vcf_file || [ -n "$p" ]
#do
#  baseInterval=$(echo ${vcf_file} | rev | cut -d'.' -f 2 | rev)
#  gatk --java-options "-Xmx128G" \
#  FilterMutectCalls \
#  -R ${genome_19} \
#  -V ${vcf_file} \
#  -O ${filter_folder}/${base}.${baseInterval}.vcf \
#  --min-median-base-quality 20 \
#  --min-median-mapping-quality 20
#done < ${interval_list}
#
## Merge all the filtered results
#count=1
#for filterFile in `ls -l ${filter_folder} | grep -E "vcf$" | awk '{print  $9}' | uniq`
#do
#  if [[ $count -eq 1 ]]
#  then 
#    cat ${filter_folder}/${filterFile} > ${filter_Mutect}
#  else
#    cat ${filter_folder}/${filterFile} | grep -vE "^#" >> ${filter_Mutect}
#  fi
#  count=`expr $count + 1`
#done
##########################################################################################

##########################################################################################
# 先整合再过滤的思路走不通，因为GATK太多bug了
#Merge output vcfs
# Error occurs 
#picard MergeVcfs \
#I=${interval_list} \
#O=${output_Mutect}

## Merge ouput vcfs by my script
#cat ${output_Mutect_gz}.${baseInterval} | grep -E "^#" > ${output_Mutect}
#while IFS="" read -r vcf_file || [ -n "$p" ]
#do
#  cat ${vcf_file} | grep -vE "^#" >> ${output_Mutect}
#done < ${interval_list}
#bgzip -f ${output_Mutect}
## Index the merged vcf
#gatk --java-options "-Xmx128G" \
#IndexFeatureFile \
#-I ${output_Mutect_gz}

#singularity exec $IMAGE_NAME gatk --java-options "-Xmx128G" \
#gatk --java-options "-Xmx128G" \
#FilterMutectCalls \
#-R ${genome_19} \
#-V ${output_Mutect_gz} \
#-O ${filter_Mutect} \
#--min-median-base-quality 20 \
#--min-median-mapping-quality 20
##########################################################################################

# Filter the Mutect2 result (followed by Cambridge's paper)
# get the PASS snv
#zcat ${filter_Mutect} | grep -vE "^#" | awk '$7~/PASS/ {print $0}' > ${filter_Mutect}.tmp
# AF <0.05 in normal,  >0.05 in tumor
#less ${filter_Mutect} | grep -E "^#" > ${filter_Mutect_Man}
#less ${filter_Mutect}.tmp | awk '{split($10,a,/:/);split($11,b,/:/);}a[3]+0 < 0.05 && b[3]+0 > 0.05{print $0}' >> ${filter_Mutect_Man}
#rm -fr ${filter_Mutect}.tmp

# Filter the Mutect2 result for the CC_data
# get the PASS snv
#cat ${filter_Mutect} | grep -vE "^#" | awk '$7~/PASS/ {print $0}' > ${filter_Mutect}.tmp
#zcat ${filter_Mutect} | grep -vE "^#" > ${filter_Mutect}.tmp
# Allele Reads in normal < 2 OR Allele Fraction in normal < 0.03
#less ${filter_Mutect}.tmp | awk '{n=split($10,a,/:/);split(a[2],b,/,/);}b[2]+0 < 2 {print $0}' > ${filter_Mutect}.tmp2
#less ${filter_Mutect}.tmp | awk '{n=split($10,a,/:/);split(a[2],b,/,/);} b[2]/a[4] < 0.03 {print $0}' > ${filter_Mutect}.tmp2
# Allele Fraction in filtered tumor reads > 0.05, the original AF in vcf is somehow wrong.
#less ${filter_Mutect}.tmp2 | awk '{n=split($11,a,/:/);split(a[2],b,/,/);} b[2]/a[4] > 0.05 {print $0}' > ${filter_Mutect}.tmp3
#less ${filter_Mutect}.tmp3 | awk '{n=split($11,a,/:/);split(a[2],b,/,/);} b[2]+0 >= 5 {print $0}' > ${filter_Mutect}.tmp4
#less ${filter_Mutect} | grep -E "^#" > ${filter_Mutect_Man}
#less ${filter_Mutect}.tmp4  >> ${filter_Mutect_Man}
#rm -fr ${filter_Mutect}.tmp*

# Get PASS result
cat ${filter_Mutect} | grep -vE "^#" | awk '$7~/PASS/ {print $0}' > ${filter_Mutect}.tmp
less ${filter_Mutect} | grep -E "^#" > ${filter_Mutect_Man}
less ${filter_Mutect}.tmp  >> ${filter_Mutect_Man}

# Annotation
/dssg/home/acct-medkwf/medkwf4/software/annovar/convert2annovar.pl --format vcf4old --includeinfo --withzyg ${filter_Mutect_Man} > ${annovar_aviinput_Mutect2}
/dssg/home/acct-medkwf/medkwf4/software/annovar/table_annovar.pl ${annovar_aviinput_Mutect2} /dssg/home/acct-medkwf/medkwf4/software/annovar/humandb -buildver hg19 -out ${annovar_annotated_Mutect2} -polish -remove -protocol refGene,icgc28,gnomad_exome,EAS.sites.2014_10,avsnp144,EAS.sites.2015_08,cosmic70,nci60  -operation g,f,f,f,f,f,f,f --otherinfo -csvout
