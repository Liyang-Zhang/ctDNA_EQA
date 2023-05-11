#!/bin/bash

#SBATCH --mail-type=end
#SBATCH --mail-user=zly18810602991@163.com
#SBATCH --job-name=test
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH -J Strelka
#SBATCH --ntasks-per-node=24
#SBATCH --output=%j.o
#SBATCH --error=%j.e

# Usage: sbatch somatic_snv_Strelka.slurm /dssg/home/acct-medkwf/medkwf4/results/MRD/CC_data/SRR3862160_CC-H007C/bam_19/SRR3862160_CC-H007C.sort.bam /dssg/home/acct-medkwf/medkwf4/results/MRD/CC_data/SRR3862162_CC-H007N/bam_19/SRR3862162_CC-H007N.sort.bam

# Initialize a variable with an intuitive name to store the name of the input fastq file

tumor=$1
normal=$2

# Get base name
baseSuper=`basename ${tumor} .sort.bam`
#base=$(echo ${baseSuper} | cut -d'_' -f 2)  # get the filename with basename function
base=${baseSuper}
echo "base is ${base}"
echo "Tumor sample is ${tumor}"
echo "Normal sample is ${normal}"

# Set working dir and other result dirs
# Remove the previous Strelka directory, otherwise it will report error
project_dir=/dssg/home/acct-medkwf/medkwf4/results/ctDNA_EQA/${base}
mkdir -p ${project_dir}

Strelka=${project_dir}/Strelka
rm -fr ${Strelka}
mkdir -p ${Strelka}

annovar=${project_dir}/ANNOVAR
mkdir -p ${annovar}

# Directory with reference genome
genome_19=/dssg/home/acct-medkwf/medkwf4/database/GATK/hg19/ucsc.hg19.fasta

# Set up environment
module load miniconda3
source activate dna

# specify the number of cores to use
cores=4

# Set output names
#Strelka_pass=${Strelka}/results/variants/${base}.vcf
annovar_aviinput_Strelka=${annovar}/${base}.Strelka
annovar_annotated_Strelka=${annovar}/${base}.Strelka.csv
Strelka_snv=${Strelka}/results/variants/somatic.snvs.vcf.gz
Strelka_indel=${Strelka}/results/variants/somatic.indels.vcf.gz
Strelka_tmp1=${Strelka}/results/variants/somatic.vcf.tmp1
Strelka_tmp2=${Strelka}/results/variants/somatic.vcf.tmp2
Strelka_filter=${Strelka}/results/variants/${base}.vcf

# Region files
# 602panel:
region_bed=/dssg/home/acct-medkwf/medkwf4/reference/bed/602panel/NanOnco_Plus_Panel_v2.0_Covered_b37_cg.parY2X.sort.chr.bed

# Use Strelka
/dssg/home/acct-medkwf/medkwf4/software/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
--normalBam ${normal} \
--tumorBam ${tumor} \
--referenceFasta ${genome_19} \
--runDir ${Strelka} \
#--exome
--callRegions ${region_bed}

${Strelka}/runWorkflow.py -m local -j 20

#zcat ${Strelka}/results/variants/somatic.snvs.vcf.gz | grep -v "^#" | awk '$7=="PASS"' > ${Strelka}/results/variants/somatic.snv.tmp.vcf
#zcat ${Strelka}/results/variants/somatic.indels.vcf.gz | grep -v "^#" | awk '$7=="PASS"' >> ${Strelka}/results/variants/somatic.snv.tmp.vcf
#zcat ${Strelka}/results/variants/somatic.snvs.vcf.gz | grep "^#" > ${Strelka_pass}
#less ${Strelka}/results/variants/somatic.snv.tmp.vcf >> ${Strelka_pass}

## Filter with the calculated Allele Fraction values
# Function making float division calculation
calc(){ awk "BEGIN { print "$*" }"; }

#split the header part in SNV file
less ${Strelka_snv} | grep -vE "^#" | awk '$7~/PASS/ {print $0}' > ${Strelka_tmp1}
less ${Strelka_snv} | grep -E "^#" > ${Strelka_filter}

# Obtain the AF information in tumor and control samples in SNV file
IFS=$'\n'
for line in `< ${Strelka_tmp1}`
do
  ref=$(echo $line | tr -d '\n' | awk '{print $4}')
  alt=$(echo $line | tr -d '\n' | awk '{print $5}')
  case $ref in

    A)
      refCountNormal=$(echo $line | tr -d '\n' | awk '{n=split($10,a,/:/);split(a[5],b,/,/);print b[1]}')
      refCounTumor=$(echo $line | tr -d '\n' | awk '{n=split($11,a,/:/);split(a[5],b,/,/);print b[1]}')
      ;;
    C)
      refCountNormal=$(echo $line | tr -d '\n' | awk '{n=split($10,a,/:/);split(a[6],b,/,/);print b[1]}')
      refCounTumor=$(echo $line | tr -d '\n' | awk '{n=split($11,a,/:/);split(a[6],b,/,/);print b[1]}')
      ;;
    G)
      refCountNormal=$(echo $line | tr -d '\n' | awk '{n=split($10,a,/:/);split(a[7],b,/,/);print b[1]}')
      refCounTumor=$(echo $line | tr -d '\n' | awk '{n=split($11,a,/:/);split(a[7],b,/,/);print b[1]}')
      ;;
    T)
      refCountNormal=$(echo $line | tr -d '\n' | awk '{n=split($10,a,/:/);split(a[8],b,/,/);print b[1]}')
      refCounTumor=$(echo $line | tr -d '\n' | awk '{n=split($11,a,/:/);split(a[8],b,/,/);print b[1]}')
      ;;
  esac

  case $alt in

    A)
      altCountNormal=$(echo $line | tr -d '\n' | awk '{n=split($10,a,/:/);split(a[5],b,/,/);print b[1]}')
      altCountTumor=$(echo $line | tr -d '\n' | awk '{n=split($11,a,/:/);split(a[5],b,/,/);print b[1]}')
      ;;
    C)
      altCountNormal=$(echo $line | tr -d '\n' | awk '{n=split($10,a,/:/);split(a[6],b,/,/);print b[1]}')
      altCountTumor=$(echo $line | tr -d '\n' | awk '{n=split($11,a,/:/);split(a[6],b,/,/);print b[1]}')
      ;;
    G)
      altCountNormal=$(echo $line | tr -d '\n' | awk '{n=split($10,a,/:/);split(a[7],b,/,/);print b[1]}')
      altCountTumor=$(echo $line | tr -d '\n' | awk '{n=split($11,a,/:/);split(a[7],b,/,/);print b[1]}')
      ;;
    T)
      altCountNormal=$(echo $line | tr -d '\n' | awk '{n=split($10,a,/:/);split(a[8],b,/,/);print b[1]}')
      altCountTumor=$(echo $line | tr -d '\n' | awk '{n=split($11,a,/:/);split(a[8],b,/,/);print b[1]}')
      ;;
  esac
  tumorDP=$((refCounTumor+altCountTumor))
  normalDP=$((refCountNormal+altCountNormal))
  tumorAF=`calc $altCountTumor/$tumorDP`
  normalAF=`calc $altCountNormal/$normalDP`
  #echo "Tumor Allele Fraction is $tumorAF, Normal Allele Fraction is $normalAF"
  if
    [ `echo "$altCountTumor > 4" | bc` -eq 1 ] && [ `echo "$tumorAF >= 0.005" | bc` -eq 1 ] && [ `echo "$normalAF < 0.03" | bc` -eq 1 ]; then
      #echo "Altered allele in tumor is $altCountTumor, tumor Allele Fraction is $tumorAF, normal Allele Fraction is $normalAF"
      echo $line >> ${Strelka_filter}
  fi
done

# filter indel result
IFS=$'\n'
less ${Strelka_indel} | grep -vE "^#" | awk '$7~/PASS/ {print $0}' > ${Strelka_tmp2}
for line in `< ${Strelka_tmp2}`
do
  refCountNormal=$(echo $line | tr -d '\n' | awk '{n=split($10,a,/:/);split(a[3],b,/,/);print b[1]}')
  altCountNormal=$(echo $line | tr -d '\n' | awk '{n=split($10,a,/:/);split(a[4],b,/,/);print b[1]}')
  refCounTumor=$(echo $line | tr -d '\n' | awk '{n=split($11,a,/:/);split(a[3],b,/,/);print b[1]}')
  altCountTumor=$(echo $line | tr -d '\n' | awk '{n=split($11,a,/:/);split(a[4],b,/,/);print b[1]}')
  tumorDP=$((refCounTumor+altCountTumor))
  normalDP=$((refCountNormal+altCountNormal))
  tumorAF=`calc $altCountTumor/$tumorDP`
  normalAF=`calc $altCountNormal/$normalDP`
  if
    [ `echo "$altCountTumor > 4" | bc` -eq 1 ] && [ `echo "$tumorAF >= 0.005" | bc` -eq 1 ] && [ `echo "$normalAF < 0.03" | bc` -eq 1 ]; then
      #echo "Altered allele in tumor is $altCountTumor, tumor Allele Fraction is $tumorAF, normal Allele Fraction is $normalAF"
      echo $line >> ${Strelka_filter}
  fi
done

# Annotation
/dssg/home/acct-medkwf/medkwf4/software/annovar/convert2annovar.pl --format vcf4old --includeinfo --withzyg ${Strelka_filter} > ${annovar_aviinput_Strelka}
/dssg/home/acct-medkwf/medkwf4/software/annovar/table_annovar.pl ${annovar_aviinput_Strelka} /dssg/home/acct-medkwf/medkwf4/software/annovar/humandb -buildver hg19 -out ${annovar_annotated_Strelka} -polish -remove -protocol refGene,icgc28,gnomad_exome,EAS.sites.2014_10,avsnp144,EAS.sites.2015_08,cosmic70,nci60  -operation g,f,f,f,f,f,f,f --otherinfo -csvout
