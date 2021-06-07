import pandas as pd
import numpy as np
import os
from inspect import getsourcefile
from os.path import abspath
import logging
logging.basicConfig(level=logging.DEBUG)

import logging
logging.basicConfig(level=logging.DEBUG)

import bin_config

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

# bwa
# def bwa(bwa, ref, fq1, fq2, )

def execmd(cmd):
    logging.INFO("cmd")
    os.system(cmd)


def seeksv(out_prefix, input_bam, ref):
    # seeksv getclip -o /path/to/outputs/prefix input.bam
    cmd1 = "{} getclip -o {} {}".format(bin_config.seeksv, out_prefix, input_bam)
    #   bwa mem  /path/to/reference.fa /path/to/prefix.clip.fq.gz | \
    #   samtools view  -Sb -o /path/to/outputs/prefix.clip.bam -
    cmd2 = "{} mem {} {}.clip.fq.gz | {} view -Sb -o {}.clip.bam".format(bin_config.bwa, out_prefix, bin_config.samtools)
    '''
    seeksv getsv /path/to/prefix.clip.bam \
             /path/to/input.bam \
             /path/to/prefix.clip.gz \
             /path/to/outputs/output.sv.txt \
             /path/to/outputs/output.unmapped.clip.fq.gz
    '''

    cmd3 = "{} getsv {}.clip.bam {} {}.clip.gz {}.sv.txt {}.unmapped.clip.fq.gz".format(bin_config.seeksv, out_prefix, input_bam, out_prefix, out_prefix, out_prefix )

    execmd(cmd1)
    execmd(cmd2)
    execmd(cmd3)
def svaba(threads, input_bam, out_dir, ref):
    # svaba run -p 64 -t srr32.sorted.markup.bam -a hpv2 -G ~/ref/hg38_hpv.fa -c 50000 -z
    prefix = os.path.join(out_dir,"svaba")
    cmd = "{} run -p {} -t {} -a {} -G {} -z".format(bin_config.svaba, threads, input_bam, prefix, ref)
    out_sv = os.path.join(out_dir,"svaba.svaba.sv.sorted.vcf.gz")
    out_txt = os.path.join(out_dir,"svaba.sv.txt")
    with open(out_txt, 'w') as f:
        sv_records = sv.read_vcf(sv_fn, precise=False)
        for i, record in enumerate(sv_records):
            f.write(str(record) + '\n')