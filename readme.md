1.generate_lh :python main.py generate_lh -f ~/remote/gate/onvirus/SRR10744032/hpv2.svaba.sv.vcf.gz.txt -b ../srr32.sorted.bam -d ~/remote/gate/onvirus/SRR10744032/srr32.depth.gz --v_chrom hpv18:1-7859 -e 3 --ploidy 2 -o test_files/svaba -s test

depth 文件： samtools depth -aa --reference ../../../gate_home/ref/hg38_hpv.fa srr32.sorted.markup.bam | bgzip -c > srr32.depth.gz && tabix -s 1 -b 2 -e 2 srr32.depth.gz

2.construct_hap： python main.py construct_hap -i test_files/svaba/test.lh -j test_files/svaba/test.junc -o test_files/svaba -s test --local_hap /home/caronkey/Documents/cityu/hpv/localHapHpv/cmake-build-debug/localHap --cbc cbc

localhap:/home/tanbowen/workspace/LocalHap/debug/src/localhap
cbc:/home/tanbowen/download/Cbc-2.9.9/build/Cbc/src/cbc

TODO seeksv 支持 done
TODO tgs 参数 done
TODO 可视化环
TODO hic筛选map	
TODO check不过

