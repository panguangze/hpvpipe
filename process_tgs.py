import sys
import os
def add_fake_lh(in_lh, out_dir):
    in_lh = open(in_lh)
    tmp_lh = os.path.join(out_dir, "tmp.lh")
    out = open(tmp_lh, "w")
    segs = []
    for line in in_lh.readlines():
        if line.startswith("SEG"):
            a = line.split(" ")
            print(a)
            segs.append(a[1].split(":")[1])
        out.write(line)
        #     in_lh.write(JUNC H:16:+ H:11:+ 969.5 -1 U B)
    for s in segs:
        out.write("JUNC H:"+s+":+ H:"+s+":+ 969.5 -1 U B\n")
        out.write("JUNC H:"+s+":+ H:"+s+":- 969.5 -1 U B\n")
        out.write("JUNC H:"+s+":- H:"+s+":- 969.5 -1 U B\n")
    in_lh.close()
    out.close()
    return tmp_lh

def parser_junc_fa_from_lh(lh_file, out_dir, ref):
    out_file = os.path.join(out_dir, "juncs.fa")
    total_len = 0
    f = open(lh_file)
    # f_o = open(out_file, "w")
    for line in f.readlines():
        a = line.split(" ")
        if a[0] == "JUNC":
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
    return out_file, total_len

