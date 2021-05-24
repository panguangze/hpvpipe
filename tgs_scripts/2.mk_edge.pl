open IN,"filter_m8.0.8";
open IN,$ARGV[0];
$len = 200;
$tr{'+'} = '-';$tr{'-'} = '+';
while(<IN>){

        chomp;
        @in = split /\t/;
        if($in[1]<=1){next}
        $infile = $in[0];

        open LS,$infile;
        %all = ();
        while(<LS>){
                chomp;
                @l = split /\t/;
                $per = $l[3]/$len;
                unless($per >0.8){next}

                @reads = split /\/|\|/,$l[1];
                $start = (split /_/,$reads[2])[0];

                @info = split /\|/,$l[0];
                $junc = $info[0];
                # $junc =~ /([H V]\d+)([+|-])\_([H V]\d+)([+|-])/;
                $junc =~ /(\d+)([+|-])\_(\d+)([+|-])/;

                ($s1,$d1,$s2,$d2) = ($1,$2,$3,$4);
                # print "$1\n";

                if($l[9] > $l[8]){
                        $minp = $start + $l[8];$dir='+';$maxp = $start + $l[9];
                }else{
                        $minp = $start + $l[9];$dir='-';$maxp = $start + $l[8];
                }
                $bp = int(($minp + $maxp)/2);

                if($dir eq '-'){
                        #if($d1 eq '+' && $d2 eq '+'){
                        #$newjunc = "$s1+ $s2+";
                        #}else{
                                $newjunc = "$s2$tr{$d2} $s1$tr{$d1}";
                                #}
                }else{
                        $newjunc = $junc;
                        $newjunc =~ s/_/ /;
                }

                $all{$bp} = "$newjunc";
        }

        #print "$infile\n";
        @junc = sort {$a<=>$b} keys  %all;
        $out = $all{$junc[0]};
        $pos = "$junc[0]";
        # print @junc;
        for($i=1;$i<@junc;$i++){
                @segs = split / /,$all{$junc[$i]};
                # print "$junc[$i]\n";
                # print "@segs\n";
                $last = (split / /,$out)[-1];
                # print "last$last\n";
                # print "out$out\n";
                if($last eq $segs[0]){
                        $out .= " $segs[1]";
                        $pos .= " $junc[$i]";
                }else{
                        print "$infile\t$out\t$pos\n";

                        $out = $all{$junc[$i]};
                        $pos = $junc[$i];
                }

        }
        print "$infile\t$out\t$pos\n";
}
