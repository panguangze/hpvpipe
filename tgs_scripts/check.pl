($segsin,$juncsin, $max_bias) = @ARGV;
open LS,$segsin;
while(<LS>){
        chomp;
        @l = split /\t/;
        $len = $l[3] - $l[2];
        $len{$l[0]} = $len;
}

open IN,$juncsin;


$tr{'-'} = '+';$tr{'+'} = '-';
while(<IN>){
        chomp;
        ($read,$juncs,$pos) = split /\t/;
        $tmp = $juncs;
        $tmp =~ s/[+|-]//g;
        @segs = split / /,$tmp;
        @pos = split / /,$pos;
        $bad = 0;
        for($i=1;$i<@pos;$i++){
                $gap = $pos[$i]-$pos[$i-1];
                $len = $len{$segs[$i]};
                #print "$_\t$segs[$i]\tGAP$gap\tLEN$len\n";
                if(abs($gap-$len) > $max_bias){
                        $bad = "$gap-$len";
                }
        }

        if($bad == 0){
                print "$juncs\n"
        }else{
                print "BAD\t$bad\t$_\n";
        }
}
