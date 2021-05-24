import pandas as pd

def process_seeksv(sv_file, threshold):
	sv = pd.read_table(sv_file, header=None,
							usecols=[0, 1, 2, 3, 4, 5, 6, 7, 9],
							names=['chrom_5p', 'pos_5p', 'strand_5p', 'left_read',
									'chrom_3p', 'pos_3p', 'strand_3p', 'right_read', 'abnormal_readpair_NO'])
	for i in range(df.shape[0]):
		if i <= df.shape[0] - 1:
			continue
		r1 = df.loc[i]
		r2 = df.loc[i+1]
		if abs(r1.start - r2.start) <= threshold and abs(r1.sv_end - r2.sv_end) <= threshold:
			l_r = max(r1.left_read, r2.left_read)
			r_r = max(r1.right_read, r2.right_read)
			ab_r = max(r1.abnormal_readpair_NO, r2.abnormal_readpair_NO)
			df.

def filter_sv_by_chrome(in_txt, out_txt, min_support_reads, chrom1, chrom2):
    in_file = open(in_txt)
    outs =[]

    for line in in_file.readlines():
        line = line.strip()
        a = line.split("\t")
        # print(a[8])
        if int(a[8]) >= min_support_read:
            # print(a[0]
			if a[0] in [chrom1, chrom2] and a[3] in [chrom1, chrom2]:
				outs.append(line)
    # print(outs)
	f_out = open(out_txt,"w")
	for l in outs:
		f_out.write(l+"\n")
	f_out.close()