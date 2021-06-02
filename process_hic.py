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

def normalize(lh_file, hic_file):
    # correct length
    h_m = pd.read_csv(hic_file,header=None, sep='\t',index_col=False, names=['s1', 's2', 'contact_v']).astype({'s1':str,'s2': str, 'contact_v': np.int64})

def seg_matrix2id_matrix(in_lh, seg_matrix_file, out_file):
    res = utils.seg2id(in_lh)
    # h_m = pd.read_csv(hic_file,header=None, sep='\t',index_col=False, names=['s1', 's2', 'contact_v']).astype({'s1':str,'s2': str, 'contact_v': np.int64})
    f_in = open(seg_matrix_file)
    # f_out = open(out_file, "w")
    for line in f_in:
        a = line.split("\t")
        id1 = res[a[0]]
        id2 = res[a[1]]
        # v = 
        # f_out.write(id1+"\t"+id2+"\t"+a[2])
    f_in.close()
    # f_out.close()
# seg_matrix2id_matrix("/home/caronkey/Documents/cityu/hpv/hpvpipe/test_files/seek_sv/seekout", "/home/caronkey/Documents/cityu/hpv/seghic.counts", "segid.counts")

# def contact2matrix()

def normalize(lh_file, seg_matrix_file):
    h_m = pd.read_csv(hic_file,header=None, sep='\t',index_col=False, names=['s1', 's2', 'contact_v']).astype({'s1':str,'s2': str, 'contact_v': np.int64})
    
    for row in h_m.itertuples():
        s1 = row["s1"].split(":")
        s1 = row["s2"].split(":")
        v = int(float(row["contact_v"]))
        l1 = int(s1[2]) - int(s1[1])
        l2 = int(s2[2]) - int(s2[1])
        v = v/(l1+l2)

