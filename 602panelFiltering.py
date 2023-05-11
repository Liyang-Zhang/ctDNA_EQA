#!/usr/bin/env python
# coding: utf-8

# In[7]:


import os
import pandas as pd
import numpy as np


# In[31]:


def get_VAF(row):
    try:
        return float(row["Otherinfo1"].split("\t")[-2].split(":")[2])
    except IndexError:
        return None

def get_TLOD(row):
    try:
        return float(row["Otherinfo1"].split("\t")[-4].split(";")[-1].split("=")[1])
    except IndexError:
        return None

def get_depth(row):
    try:
        return int(row['Otherinfo1'].split("\t")[-2].split(":")[3])
    except IndexError:
        return None

def get_read2(row):
    try:
        return int(row['Otherinfo1'].split("\t")[-2].split(":")[1].split(",")[1])
    except IndexError:
        return None

def process_file(file_path):
    data_raw = pd.read_csv(file_path)
    data_raw['VAF'] = data_raw.apply(get_VAF, axis=1)
    data_raw['TLOD'] = data_raw.apply(get_TLOD, axis=1)
    data_raw['depth'] = data_raw.apply(get_depth, axis=1)
    data_raw['mutation_reads'] = data_raw.apply(get_read2, axis=1)
    # Save the processed file with a new name
    output_file_path = file_path.replace(".csv", "_processed.csv")
    data_raw.to_csv(output_file_path, index=False)
    print(f"Processed and saved: {output_file_path}")


# In[32]:


root_dir = "/dssg/home/acct-medkwf/medkwf4/results/ctDNA_EQA"
for current_dir, subdirs, files in os.walk(root_dir):
    target_files = [file for file in files if file.endswith("hg19_multianno.csv")]
    for file in target_files:
        file_path = os.path.join(current_dir, file)
        process_file(file_path)


# In[30]:


data_raw = pd.read_csv("/dssg/home/acct-medkwf/medkwf4/results/ctDNA_EQA/202312_output-tnscope.annovar.csv.hg19_multianno.csv")
data_raw["Otherinfo1"][1].split("\t")[-2].split(":")[2]
data_raw["Otherinfo1"][1].split("\t")[-4].split(";")[-1].split("=")[1]


# In[ ]:




