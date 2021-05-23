#!/usr/bin/python
import sys
import argparse
import vcf
import re


def parser():
	parser = argparse.ArgumentParser(prog="breakpoint2vcf.py", description="Change seeksv results to vcf")
	parser.add_argument("breakpoint", help="Input breakpoint file")
	parser.add_argument("template_vcf", help="Input template vcf file")
	parser.add_argument("vcf_file", help="Output vcf file")
	return parser

base_reverse_complementary = {'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C', 'a' : 'T', 't' : 'A', 'c' : 'G', 'g' : 'C'}


def filter_bp(breakpoint_dict, breakpoint_id):	
	left_pos = int(breakpoint_dict['left_pos'])
	right_pos = int(breakpoint_dict['right_pos'])
	breakpoint_dict['junc_reads'] = 
	if breakpoint_dict['left_strand'] == '+' and breakpoint_dict['right_strand'] == '+':
		ref1 = breakpoint_dict['left_seq'][-1]
		alt1 = vcf.model._Breakend(breakpoint_dict['right_chr'], right_pos, False, True, ref1, None)
		ref2 = breakpoint_dict['right_seq'][0]
		alt2 = vcf.model._Breakend(breakpoint_dict['left_chr'], left_pos, True, False, ref2, None)
	elif breakpoint_dict['left_strand'] == '+' and breakpoint_dict['right_strand'] == '-':
		ref1 = breakpoint_dict['left_seq'][-1]
		alt1 = vcf.model._Breakend(breakpoint_dict['right_chr'], right_pos, False, False, ref1, None)
		ref2 = base_reverse_complementary[breakpoint_dict['right_seq'][0]]
		alt2 = vcf.model._Breakend(breakpoint_dict['left_chr'], left_pos, False, False, ref2, None)
	elif breakpoint_dict['left_strand'] == '-' and breakpoint_dict['right_strand'] == '+':
		ref1 = base_reverse_complementary[breakpoint_dict['left_seq'][-1]]
		alt1 = vcf.model._Breakend(breakpoint_dict['right_chr'], right_pos, True, True, ref1, None)
		ref2 = breakpoint_dict['right_seq'][0]
		alt2 = vcf.model._Breakend(breakpoint_dict['left_chr'], left_pos, True, True, ref2, None)
	
	breakend_up = "bnd" + str(breakpoint_id) + "_U"
	breakend_down = "bnd" + str(breakpoint_id) + "_D"

	record1 = vcf.model._Record(\
			CHROM=breakpoint_dict['left_chr'], \
			ID=breakend_up, \
			POS=left_pos,  \
			REF=ref1, \
			ALT=[alt1], \
			QUAL='.', \
			FILTER='PASS', \
			INFO={'SVTYPE' :'BND', 'MATEID' : breakend_down, 'CLIP_READ_NO' : breakpoint_dict['left_clip_read_NO'], \
					'STRAND' : breakpoint_dict['left_strand'], 'ABNORMAL_READPAIR_NO' : breakpoint_dict['abnormal_readpair_NO'], \
					'DEPTH' : breakpoint_dict['left_pos_depth']}, \
			FORMAT='.', \
			sample_indexes='.')

	return record1

def readBreakpointFile(breakpoint_file, template_vcf, vcf_file):
	all_records = []
	fin = open(breakpoint_file, 'r')
	header = next(fin)
	if re.match(r'@(.+)', header):
		header = header.replace('@', '')
		header = header.strip()
		header_items = header.split("\t")
	else:
		sys.stderr.write("Error: breapoint file header should start with '@'\n")
		exit(1)
	vcf_writer = vcf.Writer(open(vcf_file, 'w'), vcf.Reader(filename=template_vcf))
	breakpoint_id = 1
	for line in fin:
		line = line.strip()
		items = line.split("\t")
		breakpoint_dict = dict(zip(header_items, items))
		record1 = breakpoint2vcfRecord(breakpoint_dict, breakpoint_id)
		all_records.append(record1)
		# vcf_writer.write_record(record1)
		# vcf_writer.write_record(record2)
		breakpoint_id += 1
	
	# 
def merge_record(threshold, records):
	for i in range(len(records) - 1):
		r1 = records[i]
		r2 = records[i+1]
		if abs(r1.start - r2.start) <= threshold and abs(r1.sv_end - r2.sv_end) <= threshold:


def main():
	args = parser().parse_args(sys.argv[1:])
	readBreakpointFile(args.breakpoint, args.template_vcf, args.vcf_file)

	return 0

if __name__ == '__main__':
	main()
