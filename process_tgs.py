import sys
import os
from pyfaidx import Fasta
import collections
# add junc like 1-1, 2-2..for duplication
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
        out.write("JUNC H:"+s+":+ H:"+s+":+ 0 -1 U B\n")
        out.write("JUNC H:"+s+":+ H:"+s+":- 0 -1 U B\n")
        out.write("JUNC H:"+s+":- H:"+s+":- 0 -1 U B\n")
    in_lh.close()
    out.close()
    return tmp_lh

# 
def parser_junc_fa_from_lh(lh_file, out_dir, ref):
    ref_fa = Fasta(ref)
    out_file = os.path.join(out_dir, "juncs.fa")
    total_len = 0
    f = open(lh_file)
    # f_o = open(out_file, "w")
    id_chrom = {}
    for line in f.readlines():
        a = line.split(" ")
        if a[0] == "SEG":
            pass
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
def reverse(k):
    m = k[0:-1]
    if k[-1] == "-":
        m = m+ "+"
    if k[-1] == "+":
        m = m+ "-"
    return m

def get_segs_lens(lh_file):
    in_lh = open(lh_file)
    segs_len = {} 
    for line in in_lh.readlines():
        if line.startswith("SEG"):
            a = line.split(" ")
            infos = a[1].split(":")
            segs_len[infos[1]] = int(infos[4]) - int(infos[3])
    return segs_len
def generate_tgs_order(m8,junc_len,out_dir,lh_file, max_bias):
    segs_len = get_segs_lens(lh_file)
    tgs_out = open(os.path.join(out_dir,"tgs.juncs"),"w")
    m8_in = open(m8)
    m8_filted_r = {}
    bp_seg = {}
    final_res = {}
    for line in m8_in.readlines():
        # TODO seg lens might short than 2*len
        line_array = line.split("\t")
        if int(line_array[3]) / (2*junc_len) >= 0.8:
            tmp_junc = line_array[0].split(":")[0:2]
            # l[8] l[9]
            k = tmp_junc[0]
            v = tmp_junc[1]
            align_start = int(line_array[8])
            if int(line_array[8]) > int(line_array[9]):
                k = reverse(tmp_junc[1])
                v = reverse(tmp_junc[0])
                align_start = int(line_array[9])
            if line_array[1] in m8_filted_r.keys():
                m8_filted_r[line_array[1]].append(align_start)
            else:
                m8_filted_r[line_array[1]] = [align_start]
            bp_seg[align_start] = [k,v]
    
    for r,vs in m8_filted_r.items():
        sorted_vs = sorted(vs)
        prev_start_bp = sorted_vs[0]
        prev_segs = bp_seg[prev_start_bp]
        for i in range(1,len(vs)):
            bp = sorted_vs[i]
            segs = bp_seg[bp]
            if segs[0] == prev_segs[-1]:
                prev_segs.append(segs[1])
                if abs(segs_len[segs[0][0]] - (bp - prev_start_bp)) > segs_len[segs[0][0]] * max_bias:
                    # print(bp, prev_start_bp, segs_len[segs[0][0]] )
                    # if not 2+ 2+ duplicat, or pre_segs not greater than 2
                    if len(prev_segs) != 2 or prev_segs[0][0] == prev_segs[1][0]:
                        if " ".join(prev_segs) in final_res.keys():
                            final_res[" ".join(prev_segs)] = final_res[" ".join(prev_segs)] + 1
                        else:
                            final_res[" ".join(prev_segs)] = 1
                        # tgs_out.write(" ".join(prev_segs)+"\n")
                    prev_segs = segs
            else:
                if len(prev_segs) != 2 or prev_segs[0][0] == prev_segs[1][0]:
                    if " ".join(prev_segs) in final_res.keys():
                        final_res[" ".join(prev_segs)] = final_res[" ".join(prev_segs)] + 1
                    else:
                        final_res[" ".join(prev_segs)] = 1
                # tgs_out.write(" ".join(prev_segs)+"\n")
                prev_segs = segs
            prev_start_bp = bp
    for k,v in final_res.items():
        tgs_out.write(k+" "+str(v)+"\n")
# segs_len = get_segs_lens("/home/caronkey/Documents/cityu/hpv/hpvpipe/test_files/tgs/tmp.lh")

# if __name__ == "__main__":
#     generate_tgs_order("/home/caronkey/Documents/cityu/hpv/hpvpipe/test_files/tgs/tgs.m8", 100, "/home/caronkey/Documents/cityu/hpv/hpvpipe/test_files/tgs",segs_len, 0.15)