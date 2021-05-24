src=/home/caronkey/Documents/cityu/hpv/simulate/scripts
segs=$1
ref=$2
tgs_fa=$3
outdir=$4


cat $segs  |grep ^SEG|awk -F ':| ' '{print $3"\t"$4"\t"$5"\t"$6}' >> $outdir/input.segs
 
cat $segs |grep ^JUNC|awk -F ':| ' '{print $3" "$4" "$6" "$7}' >$outdir/input.juncs

perl $src/0.mk_junc_fa.pl $ref $outdir/input.segs $outdir/input.juncs > $outdir/juncs.fa


echo makeblastdb -in $tgs_fa -dbtype nucl -out $tgs_fa > $outdir/blast.sh
echo blastn -db $tgs_fa -query $outdir/juncs.fa -outfmt 6 \> $outdir/tgs.m8 >> $outdir/blast.sh
sh $outdir/blast.sh

perl $src/1.filter_m8.pl $outdir > $outdir/filter_m8.0.8
perl $src/2.mk_edge.pl $outdir/filter_m8.0.8  > $outdir/reads_juncs

perl $src/check.pl $outdir/input.segs $outdir/reads_juncs $seg_seq > $outdir/juncs.check
cat $outdir/juncs.check|grep -v BAD|sort|uniq > $outdir/tgs_juncs 