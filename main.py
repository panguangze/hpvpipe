import argparse
import sys, os

os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

class MainArgParser:
    def __init__(self):
        parser = argparse.ArgumentParser(prog='prelocalhap')
        parser.add_argument(dest='subfunc', help='Subcommands: ')
        args = parser.parse_args(sys.argv[1:2])
        getattr(self, args.subfunc)()

    def generate_lh(self):
        import bpsmap
        import generate_lh
        import pysam
        import pandas as pd
        import numpy as np
        parser = argparse.ArgumentParser(description='Generate localhap config for each individual')
        parser.add_argument('-f', '--sv-file',
                            dest='sv_file',
                            required=True,
                            help='Individual SV file')
        parser.add_argument('-b', '--bam-file',
                            dest='bam_file',
                            required=True,
                            help='Individual BAM file')
        parser.add_argument('-m', '--bps-map',
                            dest='bps_map',
                            required=True,
                            help='Breakpoint map file')
        parser.add_argument('-j', '--junc-db',
                            dest='junc_db',
                            required=True,
                            help='Junction database')
        parser.add_argument('-d', '--depth-tabix',
                            dest='depth_file',
                            required=True,
                            help='Tabixed depth for counting supports')
        parser.add_argument('--h_chrom',
                            dest='h_chrom',
                            required=True,
                            help='Host chrome')
        parser.add_argument('--v_chrom',
                            dest='v_chrom',
                            required=True,
                            help='Host chrome')
        parser.add_argument('-s', '--sample-name',
                            dest='sample_name',
                            required=True,
                            help='Sample name')
        parser.add_argument('-e', '--extension-bp',
                            dest='ext',
                            required=True,
                            type=int,
                            help='Extended bp for normal junctions')
        parser.add_argument('-p', '--ploidy',
                            dest='ploidy',
                            required=True,
                            default=2,
                            type=int,
                            help='Extended bp for normal junctions')
        parser.add_argument('-c', '--out-config',
                            dest='out_config',
                            required=True,
                            help='Output path of config')
        parser.add_argument('-g', '--segment',
                            dest='seg',
                            required=True,
                            help='Output path of segment')
        parser.add_argument('-i', '--keep-imprecise',
                            dest='keep_imprecise',
                            action='store_true',
                            default=True,
                            help='Keep imprecise SV')
        parser.add_argument('-I', '--keep-insertions',
                            dest='keep_insertions',
                            action='store_true',
                            default=True,
                            help='Keep insertions')
        args = parser.parse_args(sys.argv[2:])

        print('Reading SV')
        sv = bpsmap.read_sv(args.sv_file)
        sv = generate_lh.dedup(sv)

        chrome_infos = [generate_lh.parse_chrom_info(args.h_chrome), generate_lh.parse_chrom_info(args.v_chrome)]
        segs = pd.DataFrame()
        id_start = 1
        for row in chrome_infos:
            seg, id_start = generate_lh.segmentation(sv, row.chrom, row.start, row.end, id_start)
            segs = segs.append(seg)

        segs.to_csv(args.seg, index=False, sep='\t')
        bam = pysam.AlignmentFile(args.bam_file)
        depth_tabix = pysam.TabixFile(args.depth_file)

        print('Updating junc db')
        junc_db = pd.DataFrame(columns=['chrom_5p', 'pos_5p', 'strand_5p', 'chrom_3p', 'pos_3p', 'strand_3p', 'count'])
        junc_db = generate_lh.update_junc_db_by_sv(sv, junc_db)
        junc_db = generate_lh.update_junc_db_by_seg_in_chrom(segs, junc_db, bam, args.ext)
        generate_lh.write_junc_db(args.junc_db, junc_db)



        generate_lh.generate_config(args.out_config, args.sample_name, sv, segs, depth_tabix, bam, ext=args.ext, ploidy=args.ploidy, is_seeksv=args.is_seeksv)

    def parseILP(self):
        import parseILP
        parser = argparse.ArgumentParser(description='Parse ILP results')
        parser.add_argument('-i', '--in-checked-lh',
                            dest='checked_lh',
                            required=True,
                            help='Checked lh file')
        parser.add_argument('-s', '--sol',
                            dest='sol_file',
                            required=True,
                            help='Solution file')
        parser.add_argument('-o', '--output-balanced-lh',
                            dest='balanced_lh',
                            required=True,
                            help='Output balanced lh file')
        args = parser.parse_args(sys.argv[2:])

        sol = parseILP.parse_ilp_result(args.sol_file)
        parseILP.generate_balanced_lh(args.balanced_lh, args.checked_lh, sol)

    def all_process(self):
        

if __name__ == '__main__':
    MainArgParser()
