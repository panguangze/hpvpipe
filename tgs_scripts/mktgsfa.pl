$in = $ARGV[0];
$dir = $ARGV[1];
`mkdir -p $dir`;
open IN,$in;

while($head = <IN>){
        $seq = <IN>;
        chomp $head;
        @info = split />/,$head;
        open OU,">>$dir/$info[1].fa";
        print OU "$head\n$seq";

}

