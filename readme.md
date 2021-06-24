1.process_wgs: python main.py process_wgs --fq1 siha_fq1 --fq2 siha_fq2 --call_method seeksv --out_dir siha_out --ref ~/ref/hg38_hpv.fa

2.generate_lh :python main.py generate_lh -f ~/remote/hela/out.s.sv.txt -b ~/remote/hela/h8_18.bam -o test_files/hela -s hela --v_chrom HPV18REF.1:1-7857 --ploidy 3 --purity 1 --is_seeksv -e 5
3.process_tgs: python main.py process_tgs -r ~/remote/gate/ref/hg38_hpv.fa -l test_files/hela/hela.balance.lh -t /home/caronkey/Documents/cityu/hpv/tgs_fa/tgs.fa -o test_files/tgs
4.process_hic: python main.py process_hic -i test_files/hela/hela.lh --fq1 xx.f1 --fq2 xx.f2 --ref ref.fq -o test_files/hel
5.construct_hap： python main.py construct_hap -i test_files/svaba/test.lh -j test_files/svaba/test.junc -o test_files/svaba -s test --local_hap /home/caronkey/Documents/cityu/hpv/localHapHpv/cmake-build-debug/localHap --cbc cbc


