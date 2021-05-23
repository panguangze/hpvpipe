import pandas as pd

def filter_sv(sv_file, threshold):
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
		