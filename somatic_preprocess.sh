#!/bin/bash
# raw data has been instored at the following  ${project_dir} already
# USAGE : bash $0 <project directory>
# project_dir=/lustre/home/acct-medkwf/medkwf4/data/MRD/test_Pfizer
# cervical cancer dir: /lustre/home/acct-medkwf/medkwf/ncbi/public/sra/acclist/PRJNA305342_PE/raw
# ctDNA EQA tumor: /dssg/home/acct-medkwf/medkwf9/data/cgdata/DNA/cancer/ctDNA_EQA/20230420-SJZP-JZ23018508-202304231545-1
# ctDNA EQA normal: /dssg/home/acct-medkwf/medkwf9/data/cgdata/DNA/cancer/ctDNA_EQA/normal_pair
# ctDNA EQA second: /dssg/home/acct-medkwf/medkwf9/data/cgdata/DNA/cancer/ctDNA_EQA/20230505-SJZP-JZ23021426-202305092005-2

project_dir=$1

for sample_dir in `ls -l ${project_dir} | awk '{print  $9}' | grep -E "Sample_*" | uniq`
#for sample in `ls -l ${project_dir} | awk '{print  $9}' | grep -E "_1.fastq.gz$" | uniq`
do
  #echo -e "dir is ${sample_dir}"
  for sample in `ls -l ${project_dir}/${sample_dir} | awk '{print  $9}' | grep -E "_R1.fastq.gz$" | uniq`
  do
    R1=${project_dir}/${sample_dir}/${sample}
    R2=${project_dir}/${sample_dir}/${sample/_R1/_R2}
    echo -e "Sample is ${sample/_R1.fastq.gz/}"
   if
      [[ $sample == *"Lib55036"* ]]; then
      echo "running analysis on ${sample/_1.fastq.gz/}"
      sbatch /dssg/home/acct-medkwf/medkwf4/script/ctDNA_EQA/hg19/somatic_preprocess.slurm ${R1} ${R2} ${sample} ${project_dir}
   fi
  done
done

