f_d = open("/home/caronkey/Documents/cityu/hpv/hpvpipe/visual/region.tsv")
h_o = open("/home/caronkey/Downloads/Telegram Desktop/demos/forVisual2/helahpv.depth","w")
c_o = open("/home/caronkey/Downloads/Telegram Desktop/demos/forVisual2/helac.depth","w")

h_id = 1
c_id = 1
c_start = 127218091
for line in f_d.readlines():
    a = line.split("\t")
    if a[0]=="chr8":
        c_o.write("{}\t{}\t{}\t{}-{}\n".format(a[0], c_id, a[3], c_id+c_start,c_id+c_start+1))
        c_id = c_id+1
    else:
        h_o.write("{}\t{}\t{}\t{}-{}\n".format(a[0], h_id, a[3], h_id,h_id+1))
        h_id = h_id +1