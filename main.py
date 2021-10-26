import argparse
import sys, os
import utils
import bins

# os.environ["MKL_NUM_THREADS"] = "1"
# os.environ["NUMEXPR_NUM_THREADS"] = "1"
# os.environ["OMP_NUM_THREADS"] = "1"
def g_out_files(out_dir,sample,given_depth):
    suffix = ["bed","depth","lh","junc","segs"]
    res = {}
    for i in suffix:
        res[i] = os.path.join(out_dir, sample+"."+i)
    if given_depth:
        res[i] = given_depth
    return res
class MainArgParser:
    def __init__(self):
        parser = argparse.ArgumentParser(prog='prelocalhap')
        parser.add_argument(dest='subfunc', help='Subcommands: ')
        args = parser.parse_args(sys.argv[1:2])
        getattr(self, args.subfunc)()

    def generate_lh(self):
        import bpsmap
        import config
        import pysam
        import pandas as pd
        parser = argparse.ArgumentParser(description='Generate localhap config for each individual')
        parser.add_argument('--sv_file',
                            dest='sv_file',
                            required=True,
                            help='Individual SV file')
        parser.add_argument('--bam',
                            dest='bam_file',
                            required=True,
                            help='Individual BAM file')
        parser.add_argument('--v_chr',
                            dest='v_chr',
                            required=True,
                            help='virus chr name')
        parser.add_argument('--h_chrs',
                            dest='h_chrs',
                            required=True,
                            help='host chrs, if multi, split by "," eg, chr1,chr2')
        parser.add_argument('--v_len',
                            dest='v_len',
                            type=int,
                            required=True,
                            help='virus chr len')
        # parser.add_argument('--seeksv',
        #                     dest='is_seeksv',
        #                     required=False,
        #                     default=False,
        #                     action='store_true',
        #                     help='Whether seeksv results')
        parser.add_argument('--sample',
                            dest='sample',
                            required=True,
                            help='Sample name')
        parser.add_argument('--ext_bp',
                            dest='ext',
                            required=False,
                            type=int,
                            default=5,
                            help='Extended bp for normal junctions')
        parser.add_argument('--ploidy',
                            dest='ploidy',
                            required=True,
                            default=2,
                            type=int,
                            help='Extended bp for normal junctions')
        parser.add_argument('--purity',
                            dest='purity',
                            required=True,
                            default=1,
                            help='Extended bp for normal junctions')
        parser.add_argument('--out_dir',
                            dest='out_dir',
                            required=True,
                            help='Output path of config')
        parser.add_argument('--avg_whole_dp',
                            dest='avg_whole_dp',
                            required=True,
                            help='Output path of segment')
        parser.add_argument('--given_depth',
                            dest='given_depth',
                            required=False,
                            default= None,
                            help='Given depth')
        args = parser.parse_args(sys.argv[2:])
        utils.check_dir(args.out_dir)
        host_chrs = args.h_chrs.split(",")
        # print(host_chrs)
        if args.v_chr:
            host_chrs.append(args.v_chr)
        all_chrs = host_chrs
        # all output files
        out_files = g_out_files(args.out_dir,args.sample, args.given_depth)

        # generate bps map
        bps_map = bpsmap.generate_bps(args.sv_file,all_chrs,args.v_chr,args.v_len)
        # generate depth bed
        # bed_f,bam_f,depth_f,chromos,bs
        if not args.given_depth:
            bpsmap.generate_depth_bed(out_files["bed"],args.bam_file,out_files["depth"],all_chrs,bps_map)

        print('Reading SV')
        sv = pd.read_table(args.sv_file, skiprows=1, header=None,
                            usecols=[0, 1, 2, 3, 4, 5, 6, 7],
                            names=['chrom_5p', 'pos_5p', 'strand_5p', 'left_read',
                                    'chrom_3p', 'pos_3p', 'strand_3p', 'right_read'])
        config.map_bps_sv(sv, bps_map)
        # config.map_bps_chrom_infos(chrom_infos, bps_map)
        sv = config.dedup(sv)

        segs = pd.DataFrame()
        id_start = 1
        for chr in all_chrs:
            seg, id_start = config.segmentation(sv, chr,args.v_chr,args.v_len, id_start)
            segs = segs.append(seg)
        segs.to_csv(out_files["segs"], index=False, sep='\t')
        if not args.given_depth:
            bam = pysam.AlignmentFile(args.bam_file)
        depth_tabix = pysam.TabixFile(out_files["depth"])

        print('Updating junc db')
        junc_db = pd.DataFrame(columns=['chrom_5p', 'pos_5p', 'strand_5p', 'chrom_3p', 'pos_3p', 'strand_3p', 'count'])
        junc_db = config.update_junc_db_by_sv(sv, junc_db)
        junc_db = config.update_junc_db_by_seg_in_chrom(segs, junc_db, bam, args.ext)
        config.write_junc_db(out_files["junc"], junc_db)

        config.generate_config(out_files["lh"], args.sample, sv, segs, depth_tabix,bam,args.v_chr,args.avg_whole_dp, ext=args.ext,
                               ploidy=args.ploidy)

    # def generate_lh(self):
    #     import bpsmap
    #     import generate_lh
    #     import process_wgs
    #     import pysam
    #     import pandas as pd
    #     import numpy as np
    #     parser = argparse.ArgumentParser(description='Generate localhap config for each individual')
        
    #     parser.add_argument('--sv-file',
    #                         dest='sv_file',
    #                         required=True,
    #                         help='Individual SV file')
    #     parser.add_argument('--bam-file',
    #                         dest='bam_file',
    #                         required=True,
    #                         help='Individual BAM file')
    #     parser.add_argument('--out_dir',
    #                         dest='out_dir',
    #                         required=True,
    #                         help='output dir')
    #     parser.add_argument('--v_chrom',
    #                         dest='v_chrom',
    #                         required=True,
    #                         help='Virus chrome')
    #     parser.add_argument('--h_chrom',
    #                         dest='h_chrom',
    #                         required=False,
    #                         help='Host chrome, if not provide, \
    #                         we will search the main integration host chrom.')
    #     parser.add_argument('--e-bp',
    #                         dest='ext',
    #                         required=True,
    #                         type=int,
    #                         help='Extended bp for normal junctions')
    #     parser.add_argument('--ploidy',
    #                         dest='ploidy',
    #                         required=False,
    #                         default=2,
    #                         type=int,
    #                         help='Ploidy')
    #     parser.add_argument('--purity',
    #                         dest='purity',
    #                         required=False,
    #                         default=1,
    #                         type=int,
    #                         help='Ploidy')
    #     parser.add_argument('--avg_depth',
    #                         dest='avg_depth',
    #                         required=True,
    #                         help='avg_depth')
    #     parser.add_argument('--padding',
    #                         dest='padding',
    #                         required=False,
    #                         type=int,
    #                         help='padding')
    #     parser.add_argument('--is_seeksv',
    #                         dest='is_seeksv',
    #                         required=False,
    #                         action='store_true',
    #                         help='whether using span reads when calculate junction')

    #     args = parser.parse_args(sys.argv[2:])
    #     utils.check_dir(args.out_dir)
    #     host_chroms = args.h_chrs.split(",")
    #     bpsmap.generate_bps(args.sv_file,host_chroms,args.v_chr,args.v_len,args.out_dir)
    #     out_seg = os.path.join(args.out_dir, args.sample_name+'.seg')
    #     out_junc = os.path.join(args.out_dir, args.sample_name+'.junc')
    #     out_lh = os.path.join(args.out_dir, args.sample_name+'.lh')

    #     print('Reading SV')
    #     sv_sub, chrom_infos = generate_lh.filter_sv(args.sv_file,h_chrom_info, v_chrom_info, args.hic_sv, args.is_seeksv)
    #     segs = pd.DataFrame()
    #     id_start = 1
    #     for row in chrom_infos:
    #         # print(row)
    #         seg, id_start = generate_lh.segmentation(sv_sub, row['chrom'], int(row['start']), int(row['end']), id_start)
    #         segs = segs.append(seg)

    #     segs.to_csv(out_seg, index=False, sep='\t', header=False)
    #     depth_file = process_wgs.depth(args.bam_file, out_seg, args.out_dir)
    #     bam = pysam.AlignmentFile(args.bam_file)
    #     depth_tabix = pysam.TabixFile(depth_file)

    #     print('Updating junc db')
    #     junc_db = pd.DataFrame(columns=['chrom_5p', 'pos_5p', 'strand_5p', 'chrom_3p', 'pos_3p', 'strand_3p', 'count'])
    #     junc_db = generate_lh.update_junc_db_by_sv(sv_sub, junc_db, args.is_seeksv)
    #     junc_db = generate_lh.update_junc_db_by_seg_in_chrom(segs, junc_db, bam, args.ext)
    #     generate_lh.write_junc_db(out_junc, junc_db)

    #     print('Generate lh file')
    #     generate_lh.generate_config(out_lh, sv_sub, segs, depth_tabix, bam, args.is_targeted, ext=args.ext, ploidy=args.ploidy, purity=args.purity, v_chrom=v_chrom_info['chrom'],t_avg_depth=args.avg_depth, is_seeksv=args.is_seeksv)

    def generate_visual(self):
        import g_visual
        parser = argparse.ArgumentParser(description='Process tgs data')
        parser.add_argument('--hap_file',
                            dest='hap_file',
                            required=True,
                            help='Reference')
        parser.add_argument('--b_lh',
                            dest='b_lh',
                            required=True,
                            help='lh file')
        parser.add_argument('--out_dir',
                            dest='out_dir',
                            required=True,
                            help='lh file')
        g_visual.parse_hap(args.hap_file, args.b_lh, args.out_dir)
    def process_tgs(self):
        import process_tgs
        parser = argparse.ArgumentParser(description='Process tgs data')
        parser.add_argument('-r', '--ref',
                            dest='ref',
                            required=True,
                            help='Reference')
        parser.add_argument('-l','--lh_file',
                            dest='lh_file',
                            required=True,
                            help='lh file')       
        parser.add_argument('-t','--tgs_fa',
                            dest='tgs_fa',
                            required=True,
                            help='Three generation sequencing data in fasta format')
        parser.add_argument('-o','--out_dir',
                            dest='out_dir',
                            required=True,
                            help='output path for tgs file')
        parser.add_argument('--junc_len',
                            dest='junc_len',
                            default=100,
                            required=False,
                            type=int,
                            help='Length of the junction point for blast')
        parser.add_argument('--max_bias',
                            dest='max_bias',
                            default=0.15,
                            required=False,
                            type=float,
                            help='Max value of the distance between two junction point and the middle segment')
        args = parser.parse_args(sys.argv[2:])
        if not os.path.exists(args.out_dir):
            os.makedirs(args.out_dir,exist_ok=True)
        print('Parser tgs data')
        t_lh = process_tgs.add_fake_lh(args.lh_file, args.out_dir)
        tgs_cmd = "sh {}/pipe.sh {} {} {} {} {} {} {} {}".format(bins.tgs_scripts, t_lh ,args.ref,args.tgs_fa,args.out_dir,args.junc_len, args.max_bias, bins.tgs_scripts, bins.samtools)
        m8 = os.path.join(args.out_dir, "tgs.m8")
        utils.execmd(tgs_cmd)
        process_tgs.generate_tgs_order(m8,args.junc_len,args.out_dir,args.lh_file,args.max_bias)

    def process_wgs(self):
        import process_wgs
        parser = argparse.ArgumentParser(description='process wgs data')
        parser.add_argument('--fq1',
                            dest='fq1',
                            required=False,
                            help='input fq1 file')
        parser.add_argument('--fq2',
                            dest='fq2',
                            required=False,
                            help='input fq2 file')
        parser.add_argument('--ref',
                            dest='ref',
                            required=True,
                            help='input fq2 file')
        parser.add_argument('--given_bam',
                            dest='given_bam',
                            required=False,
                            help='Sorted and index bam')
        parser.add_argument('--out_dir',
                            required=True,
                            help='input fq2 file')
        parser.add_argument('--call_method',
                            required=True,
                            help='sv detection method')
        parser.add_argument('--bwa_only',
                            required=False,
                            action='store_false',
                            default=False,
                            help='Only do bwa')
        args = parser.parse_args(sys.argv[2:])
        utils.check_dir(args.out_dir)
        if args.bwa_only:
            process_wgs.bwa_wgs(args.out_dir,args.fq1, args.fq2, args.ref, args.given_bam)
            return
        if args.call_method=="seeksv":
            process_wgs.seeksv(args.out_dir, args.fq1, args.fq2, args.ref, args.given_bam)
        if args.call_method=="svaba":
            process_wgs.svaba(args.out_dir,args.fq1, args.fq2, args.ref, args.given_bam)
    def process_hic(self):
        import process_hic

        parser = argparse.ArgumentParser(description='process HIC data')
        parser.add_argument('-i', '--in_lh',
                            dest='in_lh',
                            required=True,
                            help='input lh file')
        parser.add_argument('--fq1',
                            dest='fq1',
                            required=True,
                            help='input lh file')
        parser.add_argument('--fq2',
                            dest='fq2',
                            required=True,
                            help='input lh file')
        parser.add_argument('-r', '--ref',
                            dest='ref',
                            required=True,
                            help='input junction file')
        parser.add_argument('-o', '--out_dir',
                            dest='out_dir',
                            required=True,
                            help='out_dir')
        args = parser.parse_args(sys.argv[2:])
        out_matrix = os.path.join(args.out_dir, "seg_hic.matrix")
        # Segement to fa
        seg_fa, total_len = process_hic.parser_fa_from_lh(args.in_lh, args.out_dir, args.ref)
        # bwa hic to fa
        hic_bam = process_hic.bwa_hic(args.fq1, args.fq2, seg_fa, total_len, args.out_dir)
        # bam to counts
        hic_counts = process_hic.counts(hic_bam, args.out_dir)
        # matrix normalize
        process_hic.to_matrix(args.in_lh, hic_counts, args.out_dir)
    def construct_hap(self):
        import parseILP
        parser = argparse.ArgumentParser(description='process NGS WGS data')
        parser.add_argument('-i', '--in_lh',
                            dest='in_lh',
                            required=True,
                            help='input lh file')
        parser.add_argument('-j', '--in_junc',
                            dest='in_junc',
                            required=True,
                            help='input junction file')
        parser.add_argument('-o', '--out_dir',
                            dest='out_dir',
                            required=True,
                            help='output dir')
        parser.add_argument('-s', '--sample_name',
                            dest='sample_name',
                            required=True,
                            help='sample name')
        parser.add_argument('--local_hap',
                            dest='local_hap',
                            required=True,
                            help='local hap path')
        parser.add_argument('--cbc',
                            dest='cbc_path',
                            required=True,
                            help='cbc path')
        parser.add_argument('-t','--tgs_junc',
                            dest='tgs_junc',
                            required=False,
                            help='tgs junc file')
        parser.add_argument('--is_targeted',
                            dest='is_targeted',
                            required=False,
                            default=False,
                            action='store_true',
                            help='Is your data is targeted sequencing data')

        args = parser.parse_args(sys.argv[2:])
        # call cnv
        # generate lh
        utils.check_dir(args.out_dir)
        f = os.path.join(args.out_dir, args.sample_name)
        # check
        check_cmd = "{} check {} {} {}.checked.lh {} {} --verbose".format(args.local_hap,args.in_junc, args.in_lh,f,f, args.is_targeted)
        cbc_cmd = "{} {}.lp solve solu {}.sol".format(args.cbc_path, f,f)
        utils.execmd(check_cmd)
        utils.execmd(cbc_cmd)
        sol = parseILP.parse_ilp_result(f+'.sol')
        # parser ILP result
        parseILP.generate_balanced_lh(f+'.balance.lh', f+'.checked.lh', sol)
        # generate cycle and simple haps
        solve_cmd = ""
        if args.tgs_junc:
            solve_cmd = "{} solve {} {}.balance.lh {}.circuits {}.haps --verbose".format(args.local_hap,args.in_junc,f,f,f)
        else:    
            solve_cmd = "{} solve {} {}.balance.lh {}.circuits {}.haps --verbose".format(args.local_hap,args.in_junc,f,f,f)
        utils.execmd(solve_cmd)

if __name__ == '__main__':
    MainArgParser()