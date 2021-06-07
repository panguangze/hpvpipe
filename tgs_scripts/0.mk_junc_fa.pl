# $fa = "orig_ref.fa";
# $segs = "input.segs";
($fa,$segs,$juncs,$n) = @ARGV;
open IN,$segs;
<IN>;
while(<IN>){
        chomp;
        @l = split /\t/;
        $id_seg{$l[0]} = "$l[1]:$l[2]-$l[3]";
        # print $id_seg{"H1"}
}
open IN,$juncs;
while(<IN>){
        chomp;
        ($s1,$s1d,$s2,$s2d) = split / /;

        @s1 = split /:|-/,$id_seg{$s1};
        @s2 = split /:|-/,$id_seg{$s2};

        ($p1,$p2) = ('','');
        if($s1d eq '+'){
                $t1_1 = $s1[2] - $n;
                $t1_2 = $s1[2];
        }else{
                $t1_1 = $s1[1];
                $t1_2 = $s1[1] + $n;
                $p1 .= '-i';
        }
        $fa1 = `samtools faidx $fa $p1 $s1[0]:$t1_1-$t1_2|sed 1d`;
        $fa1 =~ s/\n//g;

        if($s2d eq '+'){
                $t2_1 = $s2[1];
                $t2_2 = $s2[1] + $n;
        }else{
                $t2_1 = $s2[2] - $n;
                $t2_2 = $s2[2];
                $p2 .= '-i';
        }
        $fa2 = `samtools faidx $fa $p2 $s2[0]:$t2_1-$t2_2|sed 1d`;
        $fa2 =~ s/\n//g;

        print ">$s1$s1d\_$s2$s2d|$s1[0]:$t1_1-$t1_2\_$s2[0]:$t2_1-$t2_2\n$fa1$fa2\n";
}