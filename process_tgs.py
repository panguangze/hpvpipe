import sys
import os
def add_fake_lh(in_lh, out_dir):
    in_lh = open(in_lh)
    tmp_lh = os.path.join(out_dir, "tmp.lh")
    out = open(tmp_lh, "w")
    segs = []
    for line in in_lh.readlines():
        if line.startswith("SEG"):
            a = line.split("\t")
            segs.append(s.split(":")[1])
        out.write(line)
        #     in_lh.write(JUNC H:16:+ H:11:+ 969.5 -1 U B)
    for s in segs:
        out.write("JUNC H:"+s+":+ H:"+s+":+ 969.5 -1 U B\n")
        out.write("JUNC H:"+s+":+ H:"+s+":- 969.5 -1 U B\n")
        out.write("JUNC H:"+s+":- H:"+s+":- 969.5 -1 U B\n")
    in_lh.close()
    out.close()
    return tmp_lh