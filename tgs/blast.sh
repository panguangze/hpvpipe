makeblastdb -in /home/caronkey/remote/dl2pan/hpv_out/33ONT.fasta -dbtype nucl -out /home/caronkey/remote/dl2pan/hpv_out/33ONT.fasta
blastn -db /home/caronkey/remote/dl2pan/hpv_out/33ONT.fasta -query tgs/juncs.fa -outfmt 6 > /home/caronkey/remote/dl2pan/hpv_out/33ONT.fasta.m8
blastn -db /home/caronkey/remote/dl2pan/hpv_out/33ONT.fasta -query tgs/juncs.fa -outfmt 6 > tgs/tgs.m8
