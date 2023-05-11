# Use vcf files to  obtain mutational signatures
# author: liyang

# packages were downloaded to local conda "dna" environment
import os
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerExtractor import sigpro as sig

working_dir = "/dssg/home/acct-medkwf/medkwf4/results/MRD/CC_data/CC-H029C/mut_sig_test"
os.chdir(working_dir)

# generate matrices
#project = "mut_sig_test2"
#genome = "GRCh37"
#vcfFiles = "/dssg/home/acct-medkwf/medkwf4/results/MRD/CC_data/CC-H029C/mut_sig_test2"
#matrices = matGen.SigProfilerMatrixGeneratorFunc(project, genome, vcfFiles, exome=True, bed_file=None, chrom_based=False, plot=False, tsb_stat=False, seqInfo=False)

# Specify data type
datatype = "matrix"
data = sig.importdata(datatype)

# sigProfilerExtractor
input_type = "matrix"
output = "mutationalSignature_test"
input_data = "/dssg/home/acct-medkwf/medkwf4/results/MRD/CC_data/CC-H029C/mut_sig_test/output/SBS/test.SBS96.all"
sig.sigProfilerExtractor(input_type, output, input_data)
