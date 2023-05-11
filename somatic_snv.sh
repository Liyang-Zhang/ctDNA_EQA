#!/bin/bash
# run the variant call ready tumor and normal bqsr bam files
# USAGE : bash $0 <project directory>
# project_dir=/dssg/home/acct-medkwf/medkwf4/results/MRD/CC_data

#project_dir=$1
tool=$1
tumorID=$2
normalID=$3

#project_dir=/dssg/home/acct-medkwf/medkwf4/results/MRD/CC_data
#project_dir=/dssg/home/acct-medkwf/medkwf4/results/MRD/huashanHospital
project_dir=/dssg/home/acct-medkwf/medkwf4/results/ctDNA_EQA

for sample in `ls -l ${project_dir} | grep ^d | awk '{print  $9}' | uniq`
do
 if
  [[ $sample == ${tumorID} ]]; then
    tumor=${project_dir}/${sample}/bam_19/${sample}.sort.bam
    echo "Tumor file is ${tumor}"
 fi
 if
  [[ $sample == ${normalID} ]]; then
    normal=${project_dir}/${sample}/bam_19/${sample}.sort.bam
    echo "Normal file is ${normal}"
 fi
done

sbatch /dssg/home/acct-medkwf/medkwf4/script/ctDNA_EQA/hg19/somatic_snv_${tool}.slurm ${tumor} ${normal}

