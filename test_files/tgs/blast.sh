makeblastdb -in ../tgs_fa/33ONT.fasta -dbtype nucl -out ../tgs_fa/33ONT.fasta
blastn -db ../tgs_fa/33ONT.fasta -query test_files/tgs/juncs.fa -outfmt 6 > test_files/tgs/tgs.m8
