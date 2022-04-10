import sys
sample=sys.argv[1]
fin=sys.argv[2]

contigs=[]
count={}
c=[]
path=''
print('\t'.join(['#SampleID','haplotype_NO','colour','contig_NO','repeat_time','regid_string']))
with open(fin) as f:
        for l in f:
                if not l[0]=='#':
                        path=l.split()[-1]
print(path)
for n in path.split(','):
        if not n in c:
                c.append(n)
        else:
                c=','.join(c)
                left=c.split(n+',')[0]
                right=n+','+c.split(n+',')[-1]+','
                if not left=='':
                        contigs.append([left,'l'])
                if not right in count.keys():
                        contigs.append([right,'r'])
                        count.update({right:1})
                else:
                        count[right]+=1
                c=[n]
contigs.append([','.join(c)+',','l'])
for idx,c in enumerate(contigs):
        if c[1]=='l':
                print('\t'.join([sample,'1','#ffffff',str(idx+1),'1',c[0]]))
        else:
                print('\t'.join([sample,'1','#ffffff',str(idx+1),str(count[c[0]]),c[0]]))

 
