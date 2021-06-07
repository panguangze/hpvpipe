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
    h_m = pd.read_csv(seg_matrix_file,header=None, sep='\t',index_col=False, names=['s1', 's2', 'contact_v']).astype({'s1':str,'s2': str, 'contact_v': np.int64})
    seg_id = utils.seg2id(lh_file)
    res = []
    # print(seg_id)
    for index, row in h_m.iterrows():
        s1 = row["s1"].split(":")
        s2 = row["s2"].split(":")
        v = int(float(row["contact_v"]))
        l1 = int(s1[2]) - int(s1[1])
        l2 = int(s2[2]) - int(s2[1])
        v = v/(l1+l2)
        row["contact_v"] = v
        row["s1"] = seg_id[row["s1"]]
        row["s2"] = seg_id[row["s2"]]
        res.append(row)
    return res, len(seg_id)

def to_matrix(lh_file, seg_matrix_file, out_matrix):
    h_m, segs_len = normalize(lh_file, seg_matrix_file)
    id_copy = utils.parser_seg_info(lh_file)
    print(id_copy)
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

# def generate_matrix():
