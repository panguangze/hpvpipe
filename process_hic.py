import os
from subprocess import call
import pandas as pd
import numpy as np
import utils
import bins
import re

def parser_fa_from_lh(lh_file, out_dir, ref):
    out_file = os.path.join(out_dir, "segs.fa")
    total_len = 0
    f = open(lh_file)
    # f_o = open(out_file, "w")
    for line in f.readlines():
        a = line.split(" ")
        if a[0] == "SEG":
            b = a[1].split(":")
            s = int(b[3])
            e = int(b[4])
            total_len = total_len+(e-s)
            # m = int((s+e)/2)
            # print(b)
            cmd1 = "{} faidx {} {}:{}-{} >> {}".format(bins.samtools,ref,b[2],s,e, out_file)
            # cmd2 = "{} faidx {} {}:{}-{} >> {}".format(samtools_path,ref,b[2],m,e, out_file)
            utils.execmd(cmd1)
            # utils.execmd(cmd2)
    cmd2 = "{} index {}".format(bins.bwa, out_file)
    utils.execmd(cmd2)
    return out_file, total_len

def normalize(lh_file, hic_count_file):
    h_m = pd.read_csv(hic_count_file,header=None, sep='\t',index_col=False, names=['s1', 's2', 'contact_v']).astype({'s1':str,'s2': str, 'contact_v': np.int64})
    seg_id = utils.seg2id(lh_file)
    res = []
    # print(seg_id)
    for index, row in h_m.iterrows():
        s1 = re.split(r":|-",row["s1"])
        s2 = re.split(r":|-",row["s2"])
        # s1 = row["s1"].split(":")
        # s2 = row["s2"].split(":")
        v = int(float(row["contact_v"]))
        l1 = int(s1[2]) - int(s1[1])
        l2 = int(s2[2]) - int(s2[1])
        v = v/(l1+l2)
        row["contact_v"] = v
        row["s1"] = seg_id[row["s1"]]
        row["s2"] = seg_id[row["s2"]]
        res.append(row)
    return res, len(seg_id)

def to_matrix(lh_file, hic_count_file, out_dir):
    out_matrix = os.path.join(out_dir, "hic_matrix")
    h_m, segs_len = normalize(lh_file, hic_count_file)
    id_copy = utils.parser_seg_info(lh_file)
    # print(id_copy)
    res=np.zeros([segs_len+1,segs_len])
    for k in id_copy.keys():
        res[0][int(k)-1] = id_copy[k]
    for row in h_m:
        id1 = int(row["s1"])
        id2 = int(row["s2"])
        v = float(row["contact_v"])
        res[id1][id2-1] = v
    pd.DataFrame(res).to_csv(out_matrix, index=None, header=None)
    # return res

def bwa_hic(fq1, fq2,ref, ref_len,out_dir):
    out_bam = os.path.join(out_dir, "hic.bam")
    out_sorted_bam = os.path.join(out_dir, "hic.sorted.bam")
    # bwa mem -t 64 $1 ${var}_R1.fastq.gz ${var}_R2.fastq.gz | samtools view -@ 48 -S -h -b -F 2316 > $var.bam
    cmd1 = "{} mem -t {} {} {} {} | {} view -@ {} -S -h -b -F 2316 > {}".format(bins.bwa, bins.threads, ref, fq1, fq2, bins.samtools, bins.threads, out_bam)
    cmd2 = "{} sort -@ {} -O BAM -n -o {} {}".format(bins.samtools, bins.threads, out_sorted_bam, out_bam)
    # cmd3 = "{} index -@ {} {}".format(bins.samtools, bins.threads, out_sorted_bam)
    utils.execmd(cmd1)
    utils.execmd(cmd2)
    # utils.execmd(cmd3)
    # gc_corrected_bam = utils.gc_correction(out_sorted_bam, ref, ref_len)
    return out_sorted_bam

def counts(input_bam, out_dir):
    out_counts = os.path.join(out_dir, "hic.counts")
    # ~/apps/FALCON-Phase/bin/falcon-phase bam2 counts merged.sorted.bam merged.counts
    cmd = "{} bam2 counts {} {}".format(bins.falcon, input_bam, out_counts)
    utils.execmd(cmd)
    return out_counts
# def generate_matrix():