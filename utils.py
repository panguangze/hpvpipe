import pandas as pd
import numpy as np
import os
from inspect import getsourcefile
from os.path import abspath
import logging
logging.basicConfig(level=logging.DEBUG)

import logging
logging.basicConfig(level=logging.DEBUG)

import bins

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
    logging.INFO("cmd")
    os.system(cmd)