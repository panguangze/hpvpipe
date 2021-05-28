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
        parser.add_argument('-o', '--out_dir',
                            dest='out_dir',
                            required=True,
                            help='output dir')
        parser.add_argument('-s', '--sample_name',
                            dest='sample_name',
                            required=True,
                            help='sample_name')
        parser.add_argument('-d', '--depth-tabix',
                            dest='depth_file',
                            required=True,
                            help='Tabixed depth for counting supports')
        parser.add_argument('--v_chrom',
                            dest='v_chrom',
                            required=True,
                            help='Virus chrome')
        parser.add_argument('--h_chrom',
                            dest='h_chrom',
                            required=False,
                            help='Host chrome, if not provide, \
                            we will search the main integration host chrom.')
        parser.add_argument('-e', '--extension-bp',
                            dest='ext',
                            required=True,
                            type=int,
                            help='Extended bp for normal junctions')
        parser.add_argument('--ploidy',
                            dest='ploidy',
                            required=False,
                            default=2,
                            type=int,
                            help='Ploidy')
        parser.add_argument('--purity',
                            dest='purity',
                            required=False,
                            default=1,
                            type=int,
                            help='Ploidy')
        parser.add_argument('--is_targeted',
                            dest='is_targeted',
                            required=False,
                            default=False,
                            action='store_true',
                            help='Is your data is targeted sequencing data')
        parser.add_argument('--front_padding',
                            dest='front_padding',
                            required=False,
                            default=500,
                            type=int,
                            help='Front padding, not work if h_chrom is provided')
        parser.add_argument('--back_padding',
                            dest='back_padding',
                            required=False,
                            default=500,
                            type=int,
                            help='Back padding, not work if h_chrom is provided')
        parser.add_argument('--hic_sv',
                            dest='hic_sv',
                            required=False,
                            help='Sv file provide by hic breakfinder')


        args = parser.parse_args(sys.argv[2:])
        v_chrom_info = generate_lh.parse_chrom_info(args.v_chrom)
        h_chrom_info = ''
        if args.h_chrom:
            h_chrom_info=generate_lh.parse_chrom_info(args.h_chrom)
        else:
            h_chrom_info = generate_lh.get_integrate_chrom(args.sv_file,v_chrom_info['chrom'],args.front_padding, args.back_padding)
        out_seg = os.path.join(args.out_dir, args.sample_name+'.seg')
        out_junc = os.path.join(args.out_dir, args.sample_name+'.junc')
        out_lh = os.path.join(args.out_dir, args.sample_name+'.lh')

        print('Reading SV')
        sv_sub, chrom_infos = generate_lh.filter_sv(args.sv_file,h_chrom_info, v_chrom_info, args.hic_sv)
        segs = pd.DataFrame()
        id_start = 1
        for row in chrom_infos:
            seg, id_start = generate_lh.segmentation(sv_sub, row['chrom'], int(row['start']), int(row['end']), id_start)
            segs = segs.append(seg)

        segs.to_csv(out_seg, index=False, sep='\t')
        bam = pysam.AlignmentFile(args.bam_file)
        depth_tabix = pysam.TabixFile(args.depth_file)

        print('Updating junc db')
        junc_db = pd.DataFrame(columns=['chrom_5p', 'pos_5p', 'strand_5p', 'chrom_3p', 'pos_3p', 'strand_3p', 'count'])
        junc_db = generate_lh.update_junc_db_by_sv(sv_sub, junc_db)
        junc_db = generate_lh.update_junc_db_by_seg_in_chrom(segs, junc_db, bam, args.ext)
        generate_lh.write_junc_db(out_junc, junc_db)

        print('Generate lh file')
        generate_lh.generate_config(out_lh, sv_sub, segs, depth_tabix, bam, args.is_targeted, ext=args.ext, ploidy=args.ploidy, purity=args.purity, v_chrom=v_chrom_info['chrom'])

    def process_tgs(self):
        from subprocess import call
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
                            required=False,
                            help='Three generation sequencing data in fasta format')
        parser.add_argument('-o','--tgs_out',
                            dest='tgs_out',
                            required=True,
                            help='output path fro tgs file')
        args = parser.parse_args(sys.argv[2:])
        if not os.path.exists(args.tgs_out):
            os.mkdir(args.tgs_out)
        print('Parser tgs data')
        tgs_cmd = "sh ./tgs_scripts/pipe.sh {} {} {} {}".format(args.lh_file,args.ref,args.tgs_fa,args.tgs_out)
        if call(tgs_cmd, shell=True):
            raise Exception('Parse tgs data error')


    def construct_hap(self):
        from subprocess import call
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
        if not os.path.exists(args.out_dir):
            os.mkdir(args.out_dir)
        f = os.path.join(args.out_dir, args.sample_name)
        # check
        check_cmd = "{} check {} {} {}.checked.lh {} {} --verbose".format(args.local_hap,args.in_junc, args.in_lh,f,f, args.is_targeted)
        cbc_cmd = "{} {}.lp solve solu {}.sol".format(args.cbc_path, f,f)
        print(check_cmd)
        if call(check_cmd, shell=True):
            raise Exception("localhap check running error")
        if call(cbc_cmd, shell=True):
            raise Exception("Cbc running error")
        sol = parseILP.parse_ilp_result(f+'.sol')
        # parser ILP result
        parseILP.generate_balanced_lh(f+'.balance.lh', f+'.checked.lh', sol)
        # generate cycle and simple haps
        # TODO 调整localhap参数
        solve_cmd = ""
        if args.tgs_junc:
            solve_cmd = "{} solve {} {}.balance.lh {}.circuits {}.haps --verbose".format(args.local_hap,args.in_junc,f,f,f)
        else:    
            solve_cmd = "{} solve {} {}.balance.lh {}.circuits {}.haps --verbose".format(args.local_hap,args.in_junc,f,f,f)
        print(solve_cmd)
        if call(solve_cmd, shell=True):
            raise Exception("localhap solve running error")

if __name__ == '__main__':
    MainArgParser()