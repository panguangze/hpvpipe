import numpy as np
import pandas
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from sklearn.neighbors import KDTree
import bins
import utils


def merge_sv_tgs2sgs(sgs, tgs, thres):
    res = pd.DataFrame(columns=['chrom_5p', 'pos_5p', 'strand_5p',
                                'chrom_3p', 'pos_3p', 'strand_3p',
                                'inner_ins', 'span_reads', 'junc_reads',
                                'id', 'qual', 'filter', 'meta_info', 'anno_info'])
    res = res.astype({'chrom_5p': str, 'pos_5p': np.int64, 'strand_5p': str,
                      'chrom_3p': str, 'pos_3p': np.int64, 'strand_3p': str,
                      'inner_ins': str, 'span_reads': np.int64, 'junc_reads': np.int64,
                      'id': str, 'qual': np.float64, 'filter': str, 'meta_info': str, 'anno_info': str})
    n_found = 0
    n_found_comp = 0
    for row in tgs.itertuples():
        #         print('----------------')
        #         print(res)
        #         print('##########################')
        #         print(row)
        found = res.loc[lambda r: r.chrom_5p == row.chrom_5p] \
            .loc[lambda r: r.chrom_3p == row.chrom_3p] \
            .loc[lambda r: r.strand_5p == row.strand_5p] \
            .loc[lambda r: r.strand_3p == row.strand_3p]
        found_comp = res.loc[lambda r: r.chrom_5p == row.chrom_3p] \
            .loc[lambda r: r.chrom_3p == row.chrom_5p] \
            .loc[lambda r: r.strand_5p == ('+' if row.strand_3p == '-' else '-')] \
            .loc[lambda r: r.strand_3p == ('+' if row.strand_5p == '-' else '-')]
        #         print('+++++++++++++++++++++++++')
        #         print(found)
        if found.empty and found_comp.empty:
            found = sgs.loc[lambda r: r.chrom_5p == row.chrom_5p] \
                .loc[lambda r: r.chrom_3p == row.chrom_3p] \
                .loc[lambda r: r.strand_5p == row.strand_5p] \
                .loc[lambda r: r.strand_3p == row.strand_3p]
            found_comp = sgs.loc[lambda r: r.chrom_5p == row.chrom_3p] \
                .loc[lambda r: r.chrom_3p == row.chrom_5p] \
                .loc[lambda r: r.strand_5p == ('+' if row.strand_3p == '-' else '-')] \
                .loc[lambda r: r.strand_3p == ('+' if row.strand_5p == '-' else '-')]
        #         print('|||||||||||||||||||||||')
        #         print(found)
        #         if not found.empty or not found_comp.empty:
        #             return row, found, found_comp
        if len(found) > 0:
            n_found += 1
            dist = (found.pos_5p - row.pos_5p).abs() + (found.pos_3p - row.pos_3p).abs()
            if dist.min() <= thres * 2:
                #                 print('*************')
                #                 print(found)
                rf = found.loc[dist.idxmin()]
                idx1 = res.loc[lambda r: r.chrom_5p == rf.chrom_5p] \
                    .loc[lambda r: r.pos_5p == rf.pos_5p] \
                    .loc[lambda r: r.strand_5p == rf.strand_5p] \
                    .loc[lambda r: r.chrom_3p == rf.chrom_3p] \
                    .loc[lambda r: r.pos_3p == rf.pos_3p] \
                    .loc[lambda r: r.strand_3p == rf.strand_3p].index
                idx1_comp = res.loc[lambda r: r.chrom_5p == rf.chrom_3p] \
                    .loc[lambda r: r.pos_5p == rf.pos_3p] \
                    .loc[lambda r: r.strand_5p == ('+' if rf.strand_3p == '-' else '-')] \
                    .loc[lambda r: r.chrom_3p == rf.chrom_5p] \
                    .loc[lambda r: r.pos_3p == rf.pos_5p] \
                    .loc[lambda r: r.strand_3p == ('+' if rf.strand_5p == '-' else '-')].index
                if idx1.empty and idx1_comp.empty:
                    res = res.append(rf)
                else:
                    if not idx1.empty:
                        res.at[idx1[0], 'junc_reads'] += rf.junc_reads
                    elif not idx1_comp.empty:
                        res.at[idx1_comp[0], 'junc_reads'] += rf.junc_reads
        #                 print(f'sgs: {r.chrom_5p} {r.pos_5p} {r.strand_5p} {r.chrom_3p} {r.pos_3p} {r.strand_3p}')
        elif len(found_comp) > 0:
            n_found_comp += 1
            dist = (found_comp.pos_5p - row.pos_3p).abs() + (found_comp.pos_3p - row.pos_5p).abs()
            if dist.min() <= thres * 2:
                rf = found_comp.loc[dist.idxmin()]
                idx1 = res.loc[lambda r: r.chrom_5p == rf.chrom_5p] \
                    .loc[lambda r: r.pos_5p == rf.pos_5p] \
                    .loc[lambda r: r.strand_5p == rf.strand_5p] \
                    .loc[lambda r: r.chrom_3p == rf.chrom_3p] \
                    .loc[lambda r: r.pos_3p == rf.pos_3p] \
                    .loc[lambda r: r.strand_3p == rf.strand_3p].index
                idx1_comp = res.loc[lambda r: r.chrom_5p == rf.chrom_3p] \
                    .loc[lambda r: r.pos_5p == rf.pos_3p] \
                    .loc[lambda r: r.strand_5p == ('+' if rf.strand_3p == '-' else '-')] \
                    .loc[lambda r: r.chrom_3p == rf.chrom_5p] \
                    .loc[lambda r: r.pos_3p == rf.pos_5p] \
                    .loc[lambda r: r.strand_3p == ('+' if rf.strand_5p == '-' else '-')].index
                if idx1.empty and idx1_comp.empty:
                    res = res.append(rf)
                else:
                    if not idx1.empty:
                        res.at[idx1[0], 'junc_reads'] += rf.junc_reads
                    elif not idx1_comp.empty:
                        res.at[idx1_comp[0], 'junc_reads'] += rf.junc_reads
        else:
            rf = row
            idx1 = res.loc[lambda r: r.chrom_5p == rf.chrom_5p] \
                .loc[lambda r: r.pos_5p == rf.pos_5p] \
                .loc[lambda r: r.strand_5p == rf.strand_5p] \
                .loc[lambda r: r.chrom_3p == rf.chrom_3p] \
                .loc[lambda r: r.pos_3p == rf.pos_3p] \
                .loc[lambda r: r.strand_3p == rf.strand_3p].index
            idx1_comp = res.loc[lambda r: r.chrom_5p == rf.chrom_3p] \
                .loc[lambda r: r.pos_5p == rf.pos_3p] \
                .loc[lambda r: r.strand_5p == ('+' if rf.strand_3p == '-' else '-')] \
                .loc[lambda r: r.chrom_3p == rf.chrom_5p] \
                .loc[lambda r: r.pos_3p == rf.pos_5p] \
                .loc[lambda r: r.strand_3p == ('+' if rf.strand_5p == '-' else '-')].index
            if idx1.empty and idx1_comp.empty:
                od = rf._asdict()
                name = od['Index']
                od.pop('Index')
                ods = pd.Series(od)
                ods.name = name
                res = res.append(ods)
            else:
                if not idx1.empty:
                    res.at[idx1[0], 'junc_reads'] += rf.junc_reads
                elif not idx1_comp.empty:
                    res.at[idx1_comp[0], 'junc_reads'] += rf.junc_reads

    #                 print(f'sgs_comp: {r.chrom_5p} {r.pos_5p} {r.strand_5p} {r.chrom_3p} {r.pos_3p} {r.strand_3p}')

    return n_found, n_found_comp, res


def concat_sv(sv_list_filename):
    df = pd.DataFrame()
    with open(sv_list_filename, 'r') as fin:
        for line in fin:
            df = df.append(read_sv(line[:-1]))
    return df

def parse_sur(file_name):
    f_in = open(file_name)
    res_list = []
    for line in f_in.readlines():
        infos = line.split(" ")
        p5 = infos[1].split(":")
        p3 = infos[2].split(":")
        chrom_5p = p5[0]
        strand_5p = p5[1][0]
        pos_5p = int(p5[1][1:])

        chrom_3p = p3[0]
        strand_3p = p3[1][0]
        pos_3p = int(p3[1][1:])
        s_reads = infos[3].split("=")[1]
        res_list.append([chrom_5p, pos_5p, strand_5p, chrom_3p, pos_3p, strand_3p, int(s_reads)])
    return pandas.DataFrame(res_list, columns=['chrom_5p', 'pos_5p', 'strand_5p',
                                               'chrom_3p', 'pos_3p', 'strand_3p','junc_reads'])

def read_sv(file_name):
    return pd.read_csv(file_name, header=None, sep='\t',
                       names=['chrom_5p', 'pos_5p', 'strand_5p',
                              'chrom_3p', 'pos_3p', 'strand_3p',
                              'inner_ins', 'span_reads', 'junc_reads',
                              'id', 'qual', 'filter', 'group', 'meta_info', 'gene_5p', 'gene_3p'])


def get_precise_sv(sv_df, chrom_5p=None, start_5p=None, end_5p=None,
                   chrom_3p=None, start_3p=None, end_3p=None,
                   support_thres=5):
    # depth_tabix = pysam.TabixFile(depth_filename)
    # avg_depth = get_avg_depth(depth_tabix, chrom, start, end)
    # res_df = pd.read_table(sv_filename, header=None,
    #                       names=['chrom_5p', 'pos_5p', 'strand_5p',
    #                              'chrom_3p', 'pos_3p', 'strand_3p',
    #                              'inner_ins', 'span_reads', 'junc_reads',
    #                              'id', 'qual', 'filter', 'meta_info', 'anno_info'])
    # print(next(sv_df.itertuples()).chrom_5p, chrom)
    res_df = sv_df

    if chrom_5p:
        res_df = res_df.loc[lambda row: row.chrom_5p == chrom_5p]
        if start_5p:
            res_df = res_df.loc[lambda row: row.pos_5p >= start_5p]
        if end_5p:
            res_df = res_df.loc[lambda row: row.pos_5p <= end_5p]

    if chrom_3p:
        res_df = res_df.loc[lambda row: row.chrom_3p == chrom_3p]
        if start_3p:
            res_df = res_df.loc[lambda row: row.pos_3p >= start_3p]
        if end_3p:
            res_df = res_df.loc[lambda row: row.pos_3p <= end_3p]
    return res_df


# def get_precise_sv_seeksv(sv_filename, chrom='chr6', start=28460000, end=33500000, support_thres=5):
#     sv_df = pd.read_table(sv_filename, skiprows=1, header=None,
#                           usecols=[0, 1, 2, 3, 4, 5, 6, 7],
#                           names=['chrom_5p', 'pos_5p', 'strand_5p', 'left_read',
#                                  'chrom_3p', 'pos_3p', 'strand_3p', 'right_read'])
#     res_df = sv_df.loc[lambda row: row.chrom_5p == chrom]\
#                   .loc[lambda row: row.pos_5p >= start]\
#                   .loc[lambda row: row.pos_5p <= end]\
#                   .loc[lambda row: row.chrom_3p == chrom]\
#                   .loc[lambda row: row.pos_3p >= start]\
#                   .loc[lambda row: row.pos_3p <= end]\
#                   .loc[lambda row: row.left_read >= support_thres]\
#                   .loc[lambda row: row.right_read >= support_thres]
#     return res_df

def merge_near_pos(poses, threshold):
    r = []
    r.append(poses[0])
    for i in range(1,len(poses)):
        if poses[i] - poses[i-1] <= threshold:
            pass
        else:
            r.append(poses[i])
    return r

def get_breakpoints(sv_5p, sv_3p, is_virus):
    svs = sorted(set(sv_5p.pos_5p).union(sv_3p.pos_3p))
    print(svs)
    # if not is_virus:
    #     svs.insert(0,svs[0]-500)
    #     svs.append(svs[-1]+500)
    if is_virus:
        r = split_p_from_chr(svs)
        return sum(r,[])
    else:
        if svs[0] > 1000:
            svs.insert(0,svs[0]-1000)
        svs.append(svs[-1]+1000)
        return svs
# r = merge_near_pos(svs, 6)
    return svs


def get_breakpoints_from_list(sv_list_filename, chrom, start, end, support_thres=5):
    bps_set = set()
    n = 1
    with open(sv_list_filename, 'r') as fin:
        for line in fin:
            sv_filename, depth_filename = line[:-1].split()
            n += 1
            sv = get_precise_sv(sv_filename, depth_filename, chrom, start, end, support_thres)
            bps_set = bps_set.union(get_breakpoints(sv))
    return np.array(sorted(bps_set))



def count_neighbor(arr, r=20):
    dist = squareform(pdist(arr))
    return np.array([sum(d < r) for d in dist])


def map_bps(bps, r):
    bps_rs = bps.reshape(-1, 1)
    kdt = KDTree(bps_rs)
    ns = kdt.query_radius(bps_rs, r=r)

    inters = []
    inter = set(ns[0])
    for i in range(1, len(ns)):
        if inter.intersection(set(ns[i])):
            inter = inter.union(set(ns[i]))
        else:
            inters.append(list(inter))
            inter = set(ns[i])
    inters.append(list(inter))
    bps_map = []
    for a in inters:
        p = bps_rs[a]
        c = count_neighbor(p, r)
        pivot = p[c.argmax()][0]
        for n in p:
            bps_map.append((n[0], pivot))
    return bps_map


def split_p_from_chr(t):
    results = []
    tmp = []
    for item in t:
        if len(tmp) == 0:
            tmp.append(item)
        else:
            if item - tmp[-1] >= 100000:
                results.append(tmp)
                tmp = []
            tmp.append(item)
    results.append(tmp)
    if len(results) == 1:
        return results
    f_r = []
    for i in range(0,len(results)):
        item = results[i]
        if i == 0:
            item.append(item[-1] + 500)
        elif i == len(results)-2:
            item.append(item[0] - 500)
            # item[0] = item[0] + 500
        else:
            item.append(item[-1] + 500)
            item.append(item[0] - 500)
        f_r.append(item)
    # print(f_r)
    return f_r

def generate_depth_bed(bed_f,bam_f,depth_f,chromos,bs):
    # bps_map = [(chrom, *t) for t in bps_map]
    bed_out = open(bed_f,"w")
    for chr in chromos:
        t = bs[bs['chrom'] == chr]['after']
        tmp = split_p_from_chr(t)
        for item in tmp:
            bed_out.write("{}\t{}\t{}\n".format(chr, min(item), max(item)))
    bed_out.close()
    cmd = "{} depth -aa -b {} {} | {} -c > {} && {} -s 1 -b 2 -e 2 {}"\
        .format(bins.samtools,bed_f,bam_f,bins.bgzip,depth_f,bins.tabix,depth_f)
    utils.execmd(cmd)
    return depth_f


def generate_bps(sv_file,all_chrs,v_chr,v_len):
    sv_df = pd.read_table(sv_file, skiprows=1, header=None,
                            usecols=[0, 1, 2, 3, 4, 5, 6, 7],
                            names=['chrom_5p', 'pos_5p', 'strand_5p', 'left_read',
                                    'chrom_3p', 'pos_3p', 'strand_3p', 'right_read'])
    sv_df = sv_df.astype({
        'chrom_5p':str, 'pos_5p':np.int64, 'strand_5p':str, 'left_read':str,
        'chrom_3p':str, 'pos_3p':np.int64, 'strand_3p':str, 'right_read':str,
    })
    chroms = sorted(set(all_chrs))
    bps_map = []
    for chrom in chroms:
        sv_5p = get_precise_sv(sv_df, chrom_5p=chrom)
        sv_3p = get_precise_sv(sv_df, chrom_3p=chrom)
        bps = ""
        if v_chr and chrom == v_chr:
            bps = get_breakpoints(sv_5p, sv_3p, True) + [1, int(v_len)]
        else:
            bps = get_breakpoints(sv_5p, sv_3p, False)
        bps = np.array(sorted(set(bps)))
        bps_map.extend([(chrom, *t) for t in map_bps(bps, 10)])
    print(bps_map)
    bs = pd.DataFrame(bps_map, columns=['chrom', 'before', 'after']) \
        .sort_values(by=['chrom', 'before'])
    return bs
