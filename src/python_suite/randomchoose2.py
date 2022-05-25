import os
import sys
import time
import numpy as np
#sys.path.append("/media/sarti/data/Work/AnalyticalGenomics_Group/python_suite/")
from support import *
import sankoff
import datetime
import neighborjoining

def fast_allvsall_seqid(aln, labels, alphabet, delchars={'-', '.', 'B', 'X', 'Z'}):
#    print("START")
    ts = time.time()
    real_alphabet = sorted(list(set(alphabet) - set(delchars)))
    conv = {c : real_alphabet.index(c) for c in real_alphabet}
    delscore = 50
    delstep = 25
    alnnum = np.zeros(aln.shape)
    for i in range(alnnum.shape[0]):
        for j in range(alnnum.shape[1]):
            if aln[i,j].decode('utf-8') in delchars:
                alnnum[i,j] = delscore
                delscore += delstep
            else:
                alnnum[i,j] = conv[aln[i,j].decode('utf-8')]
#    print(alnnum)
    seqid = -np.ones((aln.shape[0], aln.shape[0]))
    for i in range(aln.shape[0]):
#        print(i)
        #print(alnnum.shape, aln[i].shape)
        diff = np.abs(alnnum - alnnum[i])
        hit = np.sum((diff == 0)*1, axis=1)
        #norm = np.sum((diff<25)*1, axis=1)
        norm = aln.shape[1]
        #seqid[i] = [hit[j]/norm[j] if norm[j] else 0 for j in range(aln.shape[0])]
        seqid[i] = [hit[j]/norm for j in range(aln.shape[0])]

    timelen = time.time() - ts
#    print(timelen)
    return seqid


# Read args
aln_fn = sys.argv[1]
list_fn = sys.argv[2]
seqs_per_bin = int(sys.argv[3])
mutations_str = sys.argv[4]
if "--generate" in sys.argv:
    generate = True
else:
    generate = False

# Fixed parameters
bin_size = 90 # Time bin in days

# Read list of labeled sequences
names = []
with open(list_fn) as lf:
    for line in lf:
        names.append(line.split('|')[3])
names_s = set(names) # Corresponding set structure

# Read time labels and divide in bins
wuhan_d = "2019-12-30"
wuhan_dt = datetime.datetime.strptime(wuhan_d, '%Y-%m-%d')
date_bins = []
this_dt = wuhan_dt
download_dt = datetime.datetime.strptime("2021-06-14", '%Y-%m-%d')
nseq_in_bins = {}
i = 0
color_list = ["rgba(255, 0, 0)", "rgba(255, 165, 0)", "rgba(255, 255, 0)", "rgba(173, 255, 47)", "rgba(0, 255, 0)", "rgba(0, 255, 175)", "rgba(0, 175, 255)"]
colors = {}
legend = {}
tmpfiles = []
while this_dt < download_dt:
    this_bin = (this_dt, this_dt + datetime.timedelta(days=bin_size))
    date_bins.append(this_bin)
    nseq_in_bins[this_bin] = 0 
    colors[this_bin] = color_list[i]
    legend[this_bin] = str(this_bin[0].date()) + " - " + str(this_bin[1].date())
    this_dt += datetime.timedelta(days=bin_size)
    if generate:
        tmpfiles.append(open("tmp_bin_{0}".format(i), "w"))
    i += 1


tmpaln_fn = "tmpaln2.afa"
if generate:
    labd = {}

    # Generate the alignments for each bin and choose an equal number of seqs from each of them
    seqs = set()
    nmaxtot = (seqs_per_bin+1)*len(nseq_in_bins)
    nseqs = 0

    # Read main aln and write bin alns (streamlined for huge files)
    # Also, add the "THIS" tag to each sequence with label belonging to input list
    with open(aln_fn) as af:
        go = False
        iline = -1
        for line in af:
            iline += 1
            if iline < 2:
                continue
            if nmaxtot == nseqs and not go:
                break
            if line.startswith(">"):
                go = False
                d = line.split('|')[2]
                if d[-2:] == '00':
                    d = d[:-1] + '1'
                if d[-5:-3] == '00':
                    d = d[:-4] + '1' + d[-3:]
                dt = datetime.datetime.strptime(d, '%Y-%m-%d')
                for ib, b in enumerate(nseq_in_bins):
                    if dt >= b[0] and dt < b[1]:
                        tmpfile = tmpfiles[ib]
                        nseq_in_bins[b] += 1
                        nseqs += 1
                        go = True
                        if line.split('|')[3] not in names_s: ### TOGLI IL NOT
                            token = " #THIS#"
                        else:
                            token = ""
                        name = line.strip() + token + '\n'
                        curr_bin = b
                        break
            elif go and line.strip() not in seqs:
                tmpfile.write(name)
                tmpfile.write(line)
                seqs.add(line.strip())
                if "#THIS#" in name:
                    token = "_THIS"
                else:
                    token = ""
                lab = name.split('|')[3]+token
                labd[name] = lab
    for tmpf in tmpfiles:
        tmpf.close()

    # Legend file (for tree visualization)
    lf = open("legend.txt", 'w')
    lf.write("TREE_COLORS\nSEPARATOR TAB\nDATA\n")

    # Main temp aln file
    taf = open(tmpaln_fn, 'w')
    with open(aln_fn) as af:
        iline = -1
        for line in af:
            iline += 1
            if iline < 2:
                taf.write(line)
                continue
            else:
                break
    for ib, b in enumerate(nseq_in_bins):
        go = False
        nrec = 0
        with open("tmp_bin_{0}".format(ib)) as tf:
            chosen_il = set()
            while nrec < seqs_per_bin:
                il = -1
                for line in tf:
                    il += 1
                    if line.startswith(">") and np.random.random() < 0.001 and nrec < seqs_per_bin and il not in chosen_il:
                        go = True
                        chosen_il.add(il)
                        taf.write(line)
                        lf.write("{0}\tbranch\t{1}\tnormal\t4\n".format(labd[line], colors[b], legend[b]))
                    else:
                        if go:
                            taf.write(line)
                            nrec += 1
                        go = False
                if go:
                    taf.write(line)
                    nrec += 1
                    go = False
                if il < seqs_per_bin:
                    print("seqs_per_bin {0} is larger than number of sequences of tmp_bin_{1} {2}".format(seqs_per_bin, ib, il))
                tf.seek(0)
    taf.close()
    lf.close()

aln, labels = parser_fasta(tmpaln_fn)
times_mx = np.zeros((len(labels), len(labels)))
dt_arr = []
for ilabel, label in enumerate(labels):
    d = label.split('|')[2]
    if d[-2:] == '00':
        d = d[:-1] + '1'
    if d[-5:-3] == '00':
        d = d[:-4] + '1' + d[-3:]
    dt = datetime.datetime.strptime(d, '%Y-%m-%d')
    dt_arr.append(dt)
mindt = min(dt_arr)
maxdt = max(dt_arr)
dtspan = maxdt - mindt
dtspan = dtspan.days
dt_arr = np.array(dt_arr)
for i in range(len(labels)):
    for j in range(len(labels)):
        ddt = dt_arr[i] - dt_arr[j]
        times_mx[i,j] = np.abs(ddt.days)/dtspan

newlabels = []
for x in labels:
    if "#THIS#" in x:
        token = "_THIS"
    else:
        token = ""
    newlabels.append(x.split('|')[3]+token)

alphabet = ['L','V','I','M','C','A','G','S','T','P','F','Y','W','E','D','N','Q','K','R','H','B','Z','X','-','.']
seqid = fast_allvsall_seqid(aln, newlabels, alphabet)
seqid = (1-seqid)
#maxseqid = np.max(seqid)
#seqid /= maxseqid
#print(maxseqid)
#print(seqid)
#dm = DistanceMatrix((1-seqid)*10, [x.split('|')[3]+"#"+x.strip()[-5:-1] for x in newlabels])
#dm = DistanceMatrix(((1-seqid)*50 + times_mx), [x.split('|')[3]+"#"+x.strip()[-5:-1] for x in newlabels])
#print(times_mx)

# NEW
nwk_t = neighborjoining.internal_nj(seqid, newlabels, {newlabels[i] : "".join([x.decode('utf-8') for x in aln[i]]) for i in range(len(aln))})

# OLD
#dm = DistanceMatrix(0*times_mx+seqid, newlabels)
#nwk_t = tree.nj(dm, result_constructor=str)

# We have to root the tree
t = Tree(nwk_t, format=1)
print(t.get_ascii(show_internal=True))
#t.set_outgroup(t.get_midpoint_outgroup())
t.set_outgroup(t.search_nodes(name=newlabels[0])[0])
nunks = 0
for node in t.traverse():
    if not node.name.strip():
        node.name = 'UNK' + str(nunks).zfill(8)
        nunks += 1

mutations = [(x[0], int(x[1:-1])-1, x[-1]) for x in sorted(mutations_str.split(','), key=lambda x: int(x[1:-1]))]
print("Calculate Sankoff")
init_d = sankoff.initialize_dictionary(alphabet)
leaves = {newlabels[i] : [(aln[i,j].decode('utf-8'), init_d[aln[i,j].decode('utf-8')]) for a,j,b in mutations] for i in range(len(newlabels))}
sankoff_score_mx = sankoff.sankoff_subs_mx(alphabet)
t, ns = sankoff.sankoff(t, sankoff_score_mx, leaves, alphabet)
t.write(outfile="tree.nhx")

# STRATEGY 1: USE SANKOFF (BUT TENDS TO GO UP ON THE TREE A BIT TOO MUCH FOR US

# STRATEGY 2: CUT AT FIRST TREE WITH ONE CHILD BEING A "THIS" LEAF, OR WITH ALL "THIS" LEAVES. THIS DOES NOT GUARANTEE TO ELIMINATE ALL "THIS", BUT WE DON'T CARE
print(t)
t2 = t.copy('deepcopy')
for node in t.traverse():
    go = True
    for l in node.iter_leaves():
        if l.name[-4:] != "THIS":
#        if sum([(mutations[i][2] == leaves[l.name][i])*1 for i in range(len(mutations))]) != len(mutations):
            go = False
            break
    for c in node.children:
        if c.is_leaf() and c.name[-4:] == "THIS":
            go = True
    if go and len(list(node.iter_leaves())) > 2:
        print(node.name, t2.search_nodes(name=node.name))
        n = t2.search_nodes(name=node.name)[0]
        nup = None
        if n.name != t2.name:
            nup = n.up
        n.detach()
        if type(nup) != type(None):
            nup.delete()
        print(node.name)
        break
print(t2)
t2, ns = sankoff.sankoff(t2, sankoff_score_mx, leaves, alphabet)

for node in t2.traverse():
    nr = 0
    for ic, c in enumerate(node.sankoff):
#        print(c, mutations[ic])
        if c == mutations[ic][2]:
            nr += 1
    node.name = str(nr) + "/" + str(len(mutations))
print(t2.get_ascii(show_internal=True))
