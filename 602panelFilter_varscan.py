#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas as pd
import numpy as np


# In[2]:


def get_mutation_info_varscan(row):
    try:
        tumorAF = float(row["Otherinfo1"].split("\t")[-1].split(":")[5].split("%")[0])
        normalAF = float(row["Otherinfo1"].split("\t")[-2].split(":")[5].split("%")[0])
        tumorDepth = float(row["Otherinfo1"].split("\t")[-1].split(":")[2])
        normalDepth = float(row["Otherinfo1"].split("\t")[-2].split(":")[2])
        tumorRefReads = float(row["Otherinfo1"].split("\t")[-1].split(":")[3])
        normalRefReads = float(row["Otherinfo1"].split("\t")[-2].split(":")[3])
        tumorMutReads = float(row["Otherinfo1"].split("\t")[-1].split(":")[4])
        normalMutReads = float(row["Otherinfo1"].split("\t")[-2].split(":")[4])
        tumorGenotype = row["Otherinfo1"].split("\t")[-1].split(":")[0]
        normalGenotype = row["Otherinfo1"].split("\t")[-2].split(":")[0]
        #info_list = [tumorAF, normalAF, tumorDepth, normalDepth, tumorRefReads, normalRefReads, tumorMutReads, normalMutReads, tumorGenotype, normalGenotype]
        #return info_list
        return (tumorAF, normalAF, tumorDepth, normalDepth, tumorRefReads, normalRefReads, tumorMutReads, normalMutReads, tumorGenotype, normalGenotype)
    except IndexError:
        return None
    
def process_file(file_path, gene_list):
    data_raw = pd.read_csv(file_path, low_memory=False)
    data_raw['rating'] = -1
    # 只保留exonic和splicing的变异
    data_sub = data_raw.loc[data_raw["Func.refGene"].str.contains("exonic|splicing")]
    # 去除synonymous SNV
    data_sub = data_sub.loc[data_raw["ExonicFunc.refGene"] != "synonymous SNV"]
    # 只保留给定的基因列表的结果
    data_sub = data_sub[data_sub["Gene.refGene"].isin(gene_list)]
    # 只保留cosmic70那列有注释的结果
    #data_sub = data_sub[data_sub["cosmic70"].notnull()]
    # 对于有CLINVAR有注释的结果的数据进行打分 (1)
    data_sub['CLINSIG'] = data_sub['CLINSIG'].fillna(".")
    data_sub.loc[data_sub['CLINSIG'].str.contains('Likely pathogenic|Pathogenic|drug response'), 'rating'] = 1
    # 进行变异注释
    data_sub["mutation_info"] = data_sub.apply(get_mutation_info_varscan, axis=1)
    data_sub[["tumorAF", "normalAF", "tumorDepth", "normalDepth", "tumorRefReads", "normalRefReads", "tumorMutReads", "normalMutReads", "tumorGenotype", "normalGenotype"]] = pd.DataFrame(data_sub["mutation_info"].tolist(), index=data_sub.index)
    output_file_path = file_path.replace(".csv.hg19_multianno.csv", "_processed.csv")
    # 直接输出大于5%的变异(rating=1)，对小于0.2%的进行打分(rating=0.5)
    data_sub['tumorAF'] = data_sub['tumorAF'].fillna('-1')
    data_sub.loc[data_sub['tumorAF'] > 5, 'rating'] = 1
    data_sub.loc[data_sub['tumorAF'] < 0.2, 'rating'] = 0.5
    # 对于变异reads数小于4的打分(rating=0.5)
    data_sub['tumorMutReads'] = data_sub['tumorMutReads'].fillna('-1')
    data_sub.loc[data_sub['tumorMutReads'] < 4, 'rating'] = 0.5
    # 输出结果
    data_sub.to_csv(output_file_path, index=False)
    print(f"Processed and saved: {output_file_path}")


# In[3]:


# 读取18 gene list
with open("/dssg/home/acct-medkwf/medkwf4/results/ctDNA_EQA/filterMerge/18gene.txt", "r") as f:
    gene_list = [line.strip() for line in f.readlines() if line.strip()]

print(gene_list)


# In[4]:


root_dir = "/dssg/home/acct-medkwf/medkwf4/results/ctDNA_EQA"
for current_dir, subdirs, files in os.walk(root_dir):
    target_files = [file for file in files if file.endswith("Varscan.csv.hg19_multianno.csv")]
    for file in target_files:
        file_path = os.path.join(current_dir, file)
        process_file(file_path, gene_list)



