import numpy as np
import pandas as pd
import sys

if len(sys.argv) < 2:
  sys.exit('usage: view depth of bam | python this.py <depth>')

depth=sys.argv[1]

depth_of_site=[]
with open(depth) as f:
        for l in f:
                depth_of_site.append(int(l.split()[1]))

if depth_of_site[0]==1:
        depth_of_site.remove(1)

sum=0
idx=0
ref_length=0
node=[]
for l in sys.stdin:
#       print str(idx)+'/'+str(len(depth_of_site))
        name=l.split()[0]
        sum=sum+int(l.split()[2])
        if idx==len(depth_of_site):
                ref_length=int(l.split()[1])
#               print ref_length
                continue
        if depth_of_site[idx]==int(l.split()[1]):
                if idx==0:
                        node.append({'name':(l.split()[0]+'_'+str(idx+1)),'weight':(sum*1.0/(depth_of_site[idx])),'site':[1,depth_of_site[idx]]})
                        sum=0
                        idx=idx+1
                else:
                        node.append({'name':(l.split()[0]+'_'+str(idx+1)),'weight':(sum*1.0/(depth_of_site[idx]+1-depth_of_site[idx-1])),'site':[depth_of_site[idx-1],depth_of_site[idx]]})
                        sum=0
                        idx=idx+1
        if idx==len(depth_of_site):
                ref_length=int(l.split()[1])
#                print ref_length

if idx==len(depth_of_site):
        node.append({'name':(l.split()[0]+'_'+str(idx+1)),'weight':(sum*1.0/(ref_length+1-depth_of_site[idx-1])),'site':[depth_of_site[idx-1],ref_length]})

#print node

for n in node:
        print n['name']+'\t'+str(n['weight'])+'\t'+str(n['site'][0])+'\t'+str(n['site'][1])+'\t'+str(n['site'][1]-n['site'][0]+1)
 
