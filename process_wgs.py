import pandas as pd
import numpy as np
import os
from inspect import getsourcefile
from os.path import abspath
import bins
from utils import execmd, gc_correction
# bwa
def bwa_wgs(out_dir,fq1, fq2, ref):
    # $bwa_bin mem -t 64 $ref $fq1 $fq2 | $samtools_bin view -bS -> $4
    out_bam = os.path.join(out_dir, "wgs.sorted.bam")
    cmd1 = "{} mem -t {} {} {} {} | {} sort -O BAM -@ {} -o {} -".format(bins.bwa, bins.threads, ref, fq1, fq2,bins.samtools,bins.threads, out_bam)
    cmd2 = "{} index {} -@ {}".format(bins.samtools, out_bam, bins.threads)
    # samtools depth -aa --reference ref_2_22_l1_real.fa L1.neo.new.bam | bgzip -c > ./L1_depth.gz && tabix -s 1 -b 2 -e 2 L1_depth.gz
    execmd(cmd1)
    execmd(cmd2)
    corrected_bam = gc_correction(out_bam, ref, bins.effectiveGenomeSize)
    return corrected_bam

def depth(input_bam, input_region, out_dir):
    out_depth=os.path.join(out_dir, "region.depth.gz")
    cmd1 = "{} depth -aa -b {} {} | {} -c > {} && {} -s 1 -b 2 -e 2 {}".format(bins.samtools, input_region, input_bam, bins.bgzip, out_depth, bins.tabix, out_depth)
    execmd(cmd1)
    return out_depth


def seeksv(out_dir, fq1, fq2, ref):
    input_bam = bwa_wgs(out_dir,fq1, fq2, ref)
    out_prefix = os.path.join(out_dir,"seeksv")
    # seeksv getclip -o /path/to/outputs/prefix input.bam
    cmd1 = "{} getclip -o {} {}".format(bins.seeksv, out_prefix, input_bam)
    #   bwa mem  /path/to/reference.fa /path/to/prefix.clip.fq.gz | \
    #   samtools view  -Sb -o /path/to/outputs/prefix.clip.bam -
    cmd2 = "{} mem {} {}.clip.fq.gz | {} view -Sb -o {}.clip.bam -".format(bins.bwa, ref, out_prefix, bins.samtools, out_prefix)
    '''
    seeksv getsv /path/to/prefix.clip.bam \
             /path/to/input.bam \
             /path/to/prefix.clip.gz \
             /path/to/outputs/output.sv.txt \
             /path/to/outputs/output.unmapped.clip.fq.gz
    '''

    cmd3 = "{} getsv {}.clip.bam {} {}.clip.gz {}.sv.txt {}.unmapped.clip.fq.gz".format(bins.seeksv, out_prefix, input_bam, out_prefix, out_prefix, out_prefix )

    execmd(cmd1)
    execmd(cmd2)
    execmd(cmd3)

def svaba(out_dir, fq1, fq2, ref):
    input_bam = bwa_wgs(out_dir,fq1, fq2, ref)
    # svaba run -p 64 -t srr32.sorted.markup.bam -a hpv2 -G ~/ref/hg38_hpv.fa -c 50000 -z
    prefix = os.path.join(out_dir,"svaba")
    cmd = "{} run -p {} -t {} -a {} -G {} -z".format(bins.svaba, bins.threads, input_bam, prefix, ref)
    execmd(cmd)
    out_sv = os.path.join(out_dir,"svaba.svaba.sv.sorted.vcf.gz")
    out_txt = os.path.join(out_dir,"svaba.sv.txt")
    with open(out_txt, 'w') as f:
        sv_records = sv.read_vcf(sv_fn, precise=False)
        for i, record in enumerate(sv_records):
            f.write(str(record) + '\n')
