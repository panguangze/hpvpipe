import pandas as pd
import numpy as np
def parser_seg_info(lh_file):
    res = []
    f_in = open("lh_file")
    for line in f_in.readlines():
        a = line.split("\t")
        if a[0] == "SEG":
            l = a[1].split(":")
            res[l[1]] = [l[2], l[3], l[4]]
    f_in.close()
    return res

def seg2id(lh_file):
    res = {}
    f_in = open(lh_file)
    for line in f_in.readlines():
        a = line.split(" ")
        if a[0] == "SEG":
            l = a[1].split(":")
            res[l[2]+":"+l[3]+"-"+l[4]] = l[1]
    f_in.close()
    print(res)
    return res