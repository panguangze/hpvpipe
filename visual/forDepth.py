import sys
import math
import numpy as np
combine_num=int(sys.argv[1])
node_margin=sys.argv[2]
orig=[]
margin=[]
for l in sys.stdin:
        orig.append(l.split())

orig=np.array(orig)
#print orig
with open(node_margin) as f:
        for l in f:
                margin.append([int(l.split()[2]),int(l.split()[3])])

seg_count=1
for m in margin:
        s=m[0]
        e=m[1]
        for i in range(int(math.ceil((e-s+1)*1.0/combine_num))):
                if (s-1+(i+1)*combine_num)<e:
                        print '\t'.join([orig[0][0],str(seg_count),str(np.mean(map(int,orig[s-1+i*combine_num:s-1+(i+1)*combine_num,2]))),str(s+i*combine_num)+'-'+str(s-1+(i+1)*combine_num)])
                else:
                        print '\t'.join([orig[0][0],str(seg_count),str(np.mean(map(int,orig[s-1+i*combine_num:e,2]))),str(s+i*combine_num)+'-'+str(e)])
                seg_count+=1
 
