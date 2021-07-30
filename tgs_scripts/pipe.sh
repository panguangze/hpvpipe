segs=$1
ref=$2
tgs_fa=$3
outdir=$4
junc_len=$5
max_bias=$6
src=$7
samtools=$8

cat $segs  |grep ^SEG|awk -F ':| ' '{print $3"\t"$4"\t"$5"\t"$6}' > $outdir/input.segs

cat $segs |grep ^JUNC|awk -F ':| ' '{print $3" "$4" "$6" "$7}' >$outdir/input.juncs

perl $src/0.mk_junc_fa.pl $ref $outdir/input.segs $outdir/input.juncs $junc_len $samtools > $outdir/juncs.fa
echo makeblastdb -in $tgs_fa -dbtype nucl -out $tgs_fa > $outdir/blast.sh
echo blastn -db $tgs_fa -query $outdir/juncs.fa -outfmt 6 \> $outdir/tgs.m8 >> $outdir/blast.sh
sh $outdir/blast.sh