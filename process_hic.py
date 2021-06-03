import os
from subprocess import call
import pandas as pd
import numpy as np
import utils
def parser_fa_from_lh(lh_file, samtools_path, out_file, ref):
    f = open(lh_file)
    # f_o = open(out_file, "w")
    for line in f.readlines():
        a = line.split(" ")
        if a[0] == "SEG":
            b = a[1].split(":")
            s = int(b[3])
            e = int(b[4])
            m = int((s+e)/2)
            # print(b)
            cmd1 = "{} faidx {} {}:{}-{} >> {}".format(samtools_path,ref,b[2],s,m, out_file)
            cmd2 = "{} faidx {} {}:{}-{} >> {}".format(samtools_path,ref,b[2],m,e, out_file)
            call(cmd1, shell=True)
            call(cmd2, shell=True)

def normalize(lh_file, seg_matrix_file):
    h_m = pd.read_csv(hic_file,header=None, sep='\t',index_col=False, names=['s1', 's2', 'contact_v']).astype({'s1':str,'s2': str, 'contact_v': np.int64})
    id_copy = utils.parser_seg_info(lh_file)
    seg_id = utils.seg2id(lh_file)
    for row in h_m.itertuples():
        s1 = row["s1"].split(":")
        s1 = row["s2"].split(":")
        copy1 = id_copy[seg_id[s1]]
        copy2 = id_copy[seg_id[s2]]
        v = int(float(row["contact_v"]))
        l1 = int(s1[2]) - int(s1[1])
        l2 = int(s2[2]) - int(s2[1])
        v = v/(l1+l2)*(copy1*copy2)
        row["contact_v"] = v
        row["s1"] = l1
        row["s2"] = l2
    return h_m, len(seg_id)

def to_matrix(lh_file, seg_matrix_file, out_matrix):
    h_m, segs_len = normalize(lh_file, seg_matrix2id_matrix)
    res=np.zeros([segs_len,segs_len])
    f_in = open(lh_file)
    for line in f_in:
        line = line.strip()
        a = line.split("\t")
        id1 = int(a[0])
        id2 = int(a[1])
        v = int(float(a[2]))
        res[id1-1][id2-1] = v
    f_in.close()
    with open('out_matrix') as f:
    for line in res:
        np.savetxt(f, line, fmt='%.2f')
    return res
