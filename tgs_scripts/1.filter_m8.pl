($indir,$len) = $ARGV;
$len = 2*$len
@file = glob("$indir/*m8");
foreach $list(@file){
        %count = ();
        $n = 0;
        open LS,$list;
        while(<LS>){
                chomp;
                @l = split /\t/;
                @info = split /-|:/,$l[0];
                #$len = $info[2]- $info[1];
                $per = $l[3]/$len;
                if($per > 0.8){
                        $count{$l[0]} ++;
                        $n ++
                }
        }
                print "$list\t$n\n";
}
