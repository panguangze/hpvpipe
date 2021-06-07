import pandas as pd
import re
import numpy as np
import bpsmap
import os
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"


def map_bps_junc(junc, bps_map):
    for i in junc.index:
        junc.at[i, 'pos_5p'] = bps_map.loc[lambda df: df.chrom == junc.loc[i, 'chrom_5p']]\
            .loc[lambda df: df.before == junc.loc[i, 'pos_5p']]\
            .iloc[0].after
        junc.at[i, 'pos_3p'] = bps_map.loc[lambda df: df.chrom == junc.loc[i, 'chrom_3p']]\
            .loc[lambda df: df.before == junc.loc[i, 'pos_3p']]\
            .iloc[0].after


def map_bps_chrom_infos(chrom_infos, bps_map):
    for i in chrom_infos.index:
        chrom_infos.at[i, 'start'] = bps_map.loc[lambda df: df.chrom == chrom_infos.loc[i, 'chrom']]\
                                            .loc[lambda df: df.before == chrom_infos.loc[i, 'start']]\
                                            .iloc[0].after
        chrom_infos.at[i, 'end'] = bps_map.loc[lambda df: df.chrom == chrom_infos.loc[i, 'chrom']]\
                                          .loc[lambda df: df.before == chrom_infos.loc[i, 'end']]\
                                          .iloc[0].after


def map_bps_sv(sv, bps_map):
    for i in sv.index:
        t = bps_map.loc[lambda df: df.chrom == sv.loc[i, 'chrom_5p']
                        ].loc[lambda df: df.before == sv.loc[i, 'pos_5p']]
        sv.at[i, 'pos_5p'] = bps_map.loc[lambda df: df.chrom == sv.loc[i, 'chrom_5p']]\
                                    .loc[lambda df: df.before == sv.loc[i, 'pos_5p']]\
                                    .iloc[0].after
        sv.at[i, 'pos_3p'] = bps_map.loc[lambda df: df.chrom == sv.loc[i, 'chrom_3p']]\
                                    .loc[lambda df: df.before == sv.loc[i, 'pos_3p']]\
                                    .iloc[0].after


def dedup(sv):
    sv = sv.sort_values(by=sv.columns[:-1].tolist())
    return sv[~sv.duplicated(sv.columns[:6], keep='last')]


def segmentation(sv, chrom, start=None, end=None, id_start=1):
    sv_5p = bpsmap.get_precise_sv(
        sv, chrom_5p=chrom, start_5p=start, end_5p=end)
    sv_3p = bpsmap.get_precise_sv(
        sv, chrom_3p=chrom, start_3p=start, end_3p=end)
    bps = sorted(set(bpsmap.get_breakpoints(sv_5p, sv_3p) + [start, end]))

    segs = []
    for p in bps[1:-1]:
        if start == None:
            start = p
            continue
        segs.append((id_start, chrom, start, p))
        start = p
        id_start += 1
    if end != None:
        segs.append((id_start, chrom, start, end))
    return pd.DataFrame(segs, columns=['ID', 'chrom', 'start', 'end']), id_start + 1


def update_junc_db_by_sv(sv, junc_db):
    for row in sv.itertuples():
        if row.inner_ins != '.':
            continue
        idx1 = junc_db.loc[lambda r: r.chrom_5p == row.chrom_5p]\
                      .loc[lambda r: r.pos_5p == row.pos_5p]\
                      .loc[lambda r: r.strand_5p == row.strand_5p]\
                      .loc[lambda r: r.chrom_3p == row.chrom_3p]\
                      .loc[lambda r: r.pos_3p == row.pos_3p]\
                      .loc[lambda r: r.strand_3p == row.strand_3p].index
        if idx1.empty:

            if True or row.junc_reads > 3:
                junc_db = junc_db.append({'chrom_5p': row.chrom_5p,
                                          'pos_5p': row.pos_5p,
                                          'strand_5p': row.strand_5p,
                                          'chrom_3p': row.chrom_3p,
                                          'pos_3p': row.pos_3p,
                                          'strand_3p': row.strand_3p,
                                          'count': 1}, ignore_index=True)
        else:
            if row.junc_reads > 5:
                junc_db.at[idx1[0], 'count'] += 1
    return junc_db


def get_normal_junc_read_num(bam, chrom, pos, ext=5):
    n = 0
    for r in bam.fetch(chrom, pos - 1, pos):
        overlapped = r.get_overlap(max(0, pos - 1 - ext), pos + ext)
        if overlapped == pos + ext - (pos - 1 - ext):
            n += 1
    return n


def update_junc_db_by_seg_in_chrom(segs, junc_db, bam, ext):
    for row in segs.iloc[:-1, :].itertuples():
        idx1 = junc_db.loc[lambda r: r.chrom_5p == row.chrom]\
                      .loc[lambda r: r.pos_5p == row.end]\
                      .loc[lambda r: r.strand_5p == '+']\
                      .loc[lambda r: r.chrom_3p == row.chrom]\
                      .loc[lambda r: r.pos_3p == row.end]\
                      .loc[lambda r: r.strand_3p == '+'].index
        if idx1.empty:
            if get_normal_junc_read_num(bam, row.chrom, row.end, ext=ext) > 3:
                junc_db = junc_db.append({'chrom_5p': row.chrom,
                                          'pos_5p': row.end,
                                          'strand_5p': '+',
                                          'chrom_3p': row.chrom,
                                          'pos_3p': row.end,
                                          'strand_3p': '+',
                                          'count': 1}, ignore_index=True)
        else:
            if get_normal_junc_read_num(bam, row.chrom, row.end, ext=ext) > 3:
                junc_db.at[idx1[0], 'count'] += 1
    return junc_db


def write_junc_db(filename, junc_db):
    junc_db.sort_values(by=['chrom_5p', 'pos_5p', 'strand_5p', 'count']).to_csv(
        filename, sep='\t', index=False)


def get_avg_depth(depth, chrom, start, end):
    return sum(map(lambda x: int(x.split('\t')[-1]), depth.fetch(chrom, start, end))) / (end - start + 1)


def generate_config(filename, sv, segs, depth_tabix, bam,is_targeted, ext, ploidy, purity, v_chrom):
    output = []
    total_depth = 0
    total_length = 0
    min_support = 3
    if is_targeted:
        min_support = 1
    with open(filename, 'w') as fout:
        output_segs = []
        for seg in segs.itertuples():
            total_length += seg.end - seg.start + 1
            seg_depth = get_avg_depth(
                depth_tabix, seg.chrom, seg.start, seg.end)
            total_depth += seg_depth * (seg.end - seg.start + 1)
            output_segs.append(f'SEG H:{seg.ID}:{seg.chrom}:{seg.start}:{seg.end} {seg_depth} -1')
        ins_id = len(segs) + 1
        ins_segs = []

        output_juncs = []
        juncs_depth = []
        left = next(segs.itertuples())
        all_supports = 0
        finded_support_num = 0
        added_right = []
        for right in segs.iloc[1:].itertuples():
            support = get_normal_junc_read_num(
                bam, left.chrom, left.end, ext=ext)
            if support > min_support:
                all_supports = all_supports + support
                finded_support_num = finded_support_num + 1
                juncs_depth.append(support)
                output_juncs.append(
                    f'JUNC H:{left.ID}:+ H:{right.ID}:+ {support} -1 U B')
                added_right.append(right.ID)
            left = right
        avg_support = all_supports / finded_support_num
        # for right in segs.iloc[1:].itertuples():
        #     print(added_right)
        #     if right.ID not in added_right:
        #         output_juncs.append(
        #             f'JUNC H:{left.ID}:+ H:{right.ID}:+ {support} -1 U B')
        #     left = right
        for row in sv.itertuples():
            if row.strand_5p == '+':
                if row.strand_3p == row.strand_5p:
                    left = segs.loc[lambda r: r.chrom == row.chrom_5p]\
                               .loc[lambda r: r.end == row.pos_5p]
                    right = segs.loc[lambda r: r.chrom == row.chrom_3p]\
                        .loc[lambda r: r.start == row.pos_3p]
                else:
                    left = segs.loc[lambda r: r.chrom == row.chrom_5p]\
                               .loc[lambda r: r.end == row.pos_5p]
                    right = segs.loc[lambda r: r.chrom == row.chrom_3p]\
                        .loc[lambda r: r.end == row.pos_3p]
            else:
                if row.strand_3p == row.strand_5p:
                    left = segs.loc[lambda r: r.chrom == row.chrom_5p]\
                               .loc[lambda r: r.start == row.pos_5p]
                    right = segs.loc[lambda r: r.chrom == row.chrom_3p]\
                        .loc[lambda r: r.end == row.pos_3p]
                else:
                    left = segs.loc[lambda r: r.chrom == row.chrom_5p]\
                               .loc[lambda r: r.start == row.pos_5p]
                    right = segs.loc[lambda r: r.chrom == row.chrom_3p]\
                        .loc[lambda r: r.start == row.pos_3p]

            juncs_depth.append(row.junc_reads)
            if row.inner_ins == '.':
                j_r = 0
                j_r = row.junc_reads
                output_juncs.append(
                    f'JUNC H:{left.ID.values[0]}:{row.strand_5p} H:{right.ID.values[0]}:{row.strand_3p} {j_r} -1 U B')

            else:
                ins_segs.append(
                    (ins_id, f'Ins_{ins_id}', 1, len(row.inner_ins), row.inner_ins))
                output_segs.append(
                    f'SEG H:{ins_id}:Ins_{ins_id}:1:{len(row.inner_ins)} 1 -1')
                output_juncs.append(
                    f'JUNC H:{left.ID.values[0]}:{row.strand_5p} H:{ins_id}:+ {row.junc_reads} -1 U B')
                output_juncs.append(
                    f'JUNC H:{ins_id}:+ H:{right.ID.values[0]}:{row.strand_3p} {row.junc_reads} -1 U B')
                ins_id += 1
        sink = ""
        for i in range(len(segs)):
            if v_chrom in segs.iloc[i].chrom.lower() and i >= 1:
                sink = segs.iloc[i-1].ID
                break
        fout.write(f'AVG_SEG_DP {total_depth * 1.0 / total_length}\n')
        fout.write(f'AVG_JUNC_DP {np.mean(juncs_depth)}\n')
        fout.write(f'PURITY {purity}\n')
        fout.write(f'AVG_PLOIDY {ploidy}\n')
        fout.write(f'PLOIDY {ploidy}m1\n')
        fout.write(f'SOURCE H:1\n')
        if is_targeted:
            fout.write(f'SINK H:{segs.iloc[-1].ID}\n')
        else:
            fout.write(f'SINK H:{sink}\n')
        fout.write('\n'.join(output_segs + output_juncs) + '\n')

    if len(ins_segs) > 0:
        pd.DataFrame(ins_segs, columns=['ID', 'chrom', 'start', 'end', 'seq']).to_csv(
            filename + '.inner_ins', index=False, sep='\t')


def parse_chrom_info(chrom_info):
    a = re.split(r'([:-])', chrom_info)
    return {'chrom': a[0], 'start': int(a[2]), 'end': int(a[4])}

def reverse_strand(strand):
    if strand == '+':
        return '-'
    if strand == '-':
        return '+'

# def filter_by_hic(sv_sub, hic_sv, ):
#     for row in sv_sub.itertuples():

def filter_by_hic(sv_sub, hic_sv, h_chrom, v_chrom):
    res = []
    hic_svs = {h_chrom:{h_chrom: [], v_chrom: []}, v_chrom: {h_chrom: [], v_chrom :[]}}
    f = open(hic_sv)
    for line in f.readlines():
        a = line.split('\t')
        chr1 = a[1]
        chr1_s = int(a[2])
        chr1_e = int(a[3])
        chr1_strand = a[4]
        chr2 = a[5]
        chr2_s = int(a[6])
        chr2_e = int(a[7])
        chr2_strand = a[8]
        
        hic_svs[chr1][chr2].append([chr1_s, chr1_e, chr1_strand, chr2_s, chr2_e, chr2_strand])
        hic_svs[chr2][chr1].append([chr2_s, chr2_e, reverse_strand(chr2_strand), chr1_s, chr1_e, reverse_strand(chr1_strand)])

    for row in sv_sub.itertuples():
        chr1 = sv_sub.chrom_5p
        chr2 = sv_sub.chrom_3p
        pos1 = sv_sub.pos_5p
        pos2 = sv_sub.pos_3p
        strand1 = sv_sub.strand_5p
        strand2 = sv_sub.strand_3p

        for l in res[chr1][chr2]:
            if l[0] <= pos1 and pos1 <= l[1] and strand1 == l[2] and l[3] <= pos2 and pos2 <= l[4] and strand2 == l[5]:
                res.append(row)
                break
            else:
                pass
    return pd.DataFrame(res)

def filter_sv(sv_file, h_chrom_info, v_chrom_info, hic_sv=None):
    res = []
    sv = bpsmap.read_sv(sv_file)
    sv = dedup(sv)
    # h_chrom_info = parse_chrom_info(h_chrom)
    # v_chrom_info = parse_chrom_info(v_chrom)
    start = int(h_chrom_info['start'])
    end = int(h_chrom_info['end'])
    # chrome_infos = [h_chrom_info,v_chrom_info]
    sv_sub = sv.loc[lambda r: (r.chrom_5p.isin([v_chrom_info['chrom'],h_chrom_info['chrom']]))
                    & (r.chrom_3p.isin([h_chrom_info['chrom'], v_chrom_info['chrom']]))]
    if hic_sv != None:
        sv_sub = filter_by_hic(sv_sub, hic_sv, h_chrom_info['chrom'], v_chrom_info['chrom'])
    for row in sv_sub.itertuples():
        if (row.chrom_5p == h_chrom_info['chrom'] \
            and (row.pos_5p > end or row.pos_5p < start))\
        or (row.chrom_3p == h_chrom_info['chrom']\
            and (row.pos_3p > end or row.pos_3p < start)):
            continue
        res.append(row)
    return pd.DataFrame(res), [h_chrom_info, v_chrom_info]

def get_integrate_chrom(sv_file, v_chrom, front_padding, back_padding):
    sv = open(sv_file)
    res={}
    for line in sv.readlines():
        a = line.split('\t')
        if a[0] == v_chrom and a[3] != v_chrom:
            tmp_p = int(a[4])
            if a[3] in res.keys():
                res[a[3]][0] = res[a[3]][0]+1
                if tmp_p < res[a[3]][1]:
                    res[a[3]][1] = tmp_p
                if tmp_p > res[a[3]][2]:
                    res[a[3]][2] = tmp_p
            else:
                res[a[3]] = []
                res[a[3]].append(1)
                res[a[3]].append(tmp_p)
                res[a[3]].append(tmp_p)
        
        if a[3] == v_chrom and a[0] != v_chrom:
            tmp_p = int(a[1])
            if a[0] in res.keys():
                res[a[0]][0] = res[a[0]][0]+1
                if tmp_p < res[a[0]][1]:
                    res[a[0]][1] = tmp_p
                if tmp_p > res[a[0]][2]:
                    res[a[0]][2] = tmp_p
            else:
                res[a[0]] = []
                res[a[0]].append(1)
                res[a[0]].append(tmp_p)
                res[a[0]].append(tmp_p)

    if len(res) == 0:
        print("There's no integration event with "+v_chrom)
        return v_chrom
    m_k = list(res.keys())[0]
    m_v = res[m_k][0]
    for k in res.keys():
        if res[k][0] >= m_v:
            m_k=k
            m_v=res[k][0]
            # TODO max染色体的长度
    s = max(res[m_k][1], 1)
    e = max(res[m_k][2], 0)
    if front_padding:
        s = max(res[m_k][1] - front_padding, 1)
    if back_padding:
        e = max(res[m_k][2] + back_padding, 0)
    return {'chrom': m_k, 'start': s, 'end': e}