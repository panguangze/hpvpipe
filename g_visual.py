import re
import os
import collections
import argparse


def split_r(hap, out_hap):
    out_f = open(out_hap, "w")
    out_f.write('\t'.join(['#Parition', 'haplotype_NO', 'colour',
                'contig_NO', 'repeat_time', 'regid_string'])+"\n")
    for k, v in hap.items():
        for h_num, path in enumerate(v):
            # print(path)
            contigs = []
            count = {}
            c = []
            # print(path)
            for n in path[1:]:
                # print(n)
                if not n in c:
                    c.append(n)
                else:
                    c = ','.join(c)
                    left = c.split(n+',')[0]
                    right = n+','+c.split(n+',')[-1]+','
                    if not left == '':
                        contigs.append([left, 'l'])
                    if not right in count.keys():
                        contigs.append([right, 'r'])
                        count.update({right: 1})
                    else:
                        count[right] += 1
                    c = [n]
            contigs.append([','.join(c)+',', 'l'])
            for idx, c in enumerate(contigs):
                print(k)
                if c[1] == 'l':
                    out_f.write('\t'.join([k, ",".join([str(i) for i in range(
                        0, int(path[0]))]), '#ffffff', str(idx+1), '1', c[0]])+'\n')
                else:
                    out_f.write('\t'.join([k, ",".join([str(i) for i in range(
                        0, int(path[0]))]), '#ffffff', str(idx+1), str(count[c[0]]), c[0]])+'\n')


def parse_balance_lh(balanced_lh, out_seg):
    out_fin = open(out_seg, "w")
    out_fin.write('chrom\tstart\tend\tid\tcopy_origin\tcopy_balanced\n')
    b_in = open(balanced_lh)
    for line in b_in.readlines():
        if "SEG " in line:
            line = line.strip()
            line = line.split(" ")
            a = line[1].split(":")
            out_fin.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(
                a[2], a[3], a[4], a[1], line[3], line[-1]))


def parse_hap(hap_file, balanced_lh, out_dir):
    hap_in = open(hap_file)
    res = {}
    for line in hap_in.readlines():
        line = line.strip("\n")
        info = line.split(":")
        hap_count = collections.Counter(info[1:])
        for hap, count in hap_count.items():
            vs = []
            vs.append(str(count))
            for i in hap.split(" ")[:-1]:
                if "-" in i:
                    vs.append(i[0:-1]+"_r")
                else:
                    vs.append(i[0:-1])
            if info[0] not in res.keys():
                res[info[0]] = [vs]
            else:
                res[info[0]].append(vs)
    print(res)
    out_seg = os.path.join(out_dir, "visual.seg")
    out_hap = os.path.join(out_dir, "visual.hap")
    split_r(res, out_hap)
    parse_balance_lh(balanced_lh, out_seg)


def main():
    parser = argparse.ArgumentParser(
        description='Generate localhap config for each individual')
    parser.add_argument('--hap',
                        dest='hap_file',
                        required=True,
                        help='haplotype file')
    parser.add_argument('--balanced_lh',
                        dest='balanced_lh',
                        required=True,
                        help='balanced_lh file')
    parser.add_argument('--out_file',
                        dest='out_file',
                        required=True,
                        help='out file')
    args = parser.parse_args()
    parse_hap(args.hap_file, args.balanced_lh, args.out_file)

if __name__ == "__main__":
    main()
