import os
from subprocess import call
def parser_fa_from_lh(lh_file, samtools_path, out_file, ref):
    f = open(lh_file)
    # f_o = open(out_file, "w")
    for line in f.readlines():
        a = line.split(" ")
        if a[0] == "SEG":
            b = a[1].split(":")
            print(b)
            cmd = "{} faidx {} {}:{}-{} >> {}".format(samtools_path,ref,b[2],b[3],b[4], out_file)
            call(cmd, shell=True)

def normalize(lh_file, hic_matrix):
    
