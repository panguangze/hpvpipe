import pandas as pd
import numpy as np
import os
from inspect import getsourcefile
from os.path import abspath
import logging
import bins
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s')
# process lh_file
def parser_seg_info(lh_file):
    res = {}
    f_in = open(lh_file)
    for line in f_in.readlines():
        a = line.split(" ")
        if a[0] == "SEG":
            l = a[1].split(":")
            res[l[1]] = a[3]
    f_in.close()
    return res

def seg2id(lh_file):
    res = {}
    f_in = open(lh_file)
    for line in f_in.readlines():
        a = line.split(" ")
        if a[0] == "SEG":
            l = a[1].split(":")
            res[l[2]+":"+l[3]+":"+l[4]] = l[1]
    f_in.close()
    # print(res)
    return res

def execmd(cmd):
    logging.info(cmd)
    # os.system(cmd)

# gc corrected
def gc_correction(input_bam, ref, effectiveGenomeSize):
    corrected_bam = input_bam.replace(".bam", ".gccorrected.bam")
    # faToTwoBit hg38_hpv.fa hg38_hpv.bit
    cmd1 = "{} {} {}.2bit".format(bins.faToTwoBit, ref, ref)
    #  computeGCBias -b file.bam --effectiveGenomeSize 2150570000 -g mm9.2bit -l 200 --GCbiasFrequenciesFile freq.txt [options]
    cmd2 = "{} -b {} --effectiveGenomeSize {} -g {}.2bit --GCbiasFrequenciesFile {}.freq.txt".format(bins.computeGCBias, input_bam, effectiveGenomeSize, ref, ref)
    #  correctGCBias -b file.bam --effectiveGenomeSize 2150570000 -g mm9.2bit --GCbiasFrequenciesFile freq.txt -o gc_corrected.bam [options]
    cmd3 = "{} -b {} --effectiveGenomeSize {} -g {}.2bit --GCbiasFrequenciesFile {}.freq.txt -o {}".format(bins.correctGCBias, input_bam, effectiveGenomeSize, ref, ref, corrected_bam)
    cmd4 = "{} index {} -@ {}".format(bins.samtools, corrected_bam, bins.threads)
    execmd(cmd1)
    execmd(cmd2)
    execmd(cmd3)
    execmd(cmd4)
    return corrected_bam