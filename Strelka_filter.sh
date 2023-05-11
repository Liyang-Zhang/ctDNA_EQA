#!/bin/bash

# Function making float division calculation
calc(){ awk "BEGIN { print "$*" }"; }

# Set file names
VCF=/lustre/home/acct-medkwf/medkwf4/results/MRD/CC_data/CC-H007C/Strelka/results/variants/somatic.snvs.vcf.gz
VCF_tmp1=/lustre/home/acct-medkwf/medkwf4/results/MRD/CC_data/CC-H007C/Strelka/results/variants/somatic.snvs.vcf.tmp1
Strelka_filter=/lustre/home/acct-medkwf/medkwf4/results/MRD/CC_data/CC-H007C/Strelka/results/variants/somatic.snvs.filter.vcf

# split the header part
less ${VCF} | grep -vE "^#" | awk '$7~/PASS/ {print $0}' > ${VCF_tmp1}
less ${VCF} | grep -E "^#" > ${Strelka_filter}

# Obtain the AF information in tumor and control samples
while read line || [[ -n $line ]];
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
    [ `echo "$altCountTumor > 4" | bc` -eq 1 ] && [ `echo "$tumorAF >= 0.05" | bc` -eq 1 ] && [ `echo "$normalAF < 0.03" | bc` -eq 1 ]; then
      #echo "Altered allele in tumor is $altCountTumor, tumor Allele Fraction is $tumorAF, normal Allele Fraction is $normalAF"
      echo $line | awk '{print $0}' >> ${Strelka_filter}
  fi
done < ${VCF_tmp1}

# similar script for indel filter

