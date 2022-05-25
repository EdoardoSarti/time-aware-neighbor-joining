import os
import sys
import time
import numpy as np
sys.path.append("/media/sarti/data/Work/AnalyticalGenomics_Group/python_suite/")
from support import *
import sankoff
import datetime
import matplotlib.pyplot as plt

"""
def errprint(text):
    print(text, file=sys.stderr)


def parser_fasta(aln_fn):
    if not os.path.exists(aln_fn):
        errprint("[parser_fasta] FATAL: input file does not exist")
        errprint("aln_fn:", aln_fn)
        exit(1)
    firstline = True
    aln_seqs = {}
    names = []
    with open(aln_fn) as af:
        for line in af:
            if firstline and not line.startswith(">"):
                errprint("[parser_fasta] FATAL: wrong format")
                errprint("line 1:", line)
                exit(1)
            if line.startswith(">"):
                name = line[1:].strip()
                if name in names:
                    errprint("[parser_fasta] FATAL: two sequences share the same header")
                    errprint("Header: {0}\nSequence 1: {1}\nSequence 2: {2}".format(name, names.index(name)+1, len(names)+1))
                    exit(1)
                names.append(name)
                aln_seqs[name] = ""
            else:
                aln_seqs[name] += line.strip()
            if firstline:
                firstline = False

    if not names:
        errprint("[parser_fasta] WARNING: empty alignment")
        return []

    L = len(aln_seqs[names[0]])
    M = len(names)
    out_aln = np.chararray((M, L))
    errprint("[parser_fasta] REMARK: filename {0}, length {1}, depth {2}".format(os.path.basename(aln_fn), L, M))

    for iname, name in enumerate(names):
        if len(aln_seqs[name]) != L:
            errprint("[parser_fasta] FATAL: sequence {0} is long {1} instead of {2}".format(iname+1, len(len(aln_seqs[name]), L)))
            exit(1)
        out_aln[iname] = list(aln_seqs[name])
    return out_aln, names
"""

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
    print(alnnum)
    seqid = -np.ones((aln.shape[0], aln.shape[0]))
    for i in range(aln.shape[0]):
#        print(i)
        #print(alnnum.shape, aln[i].shape)
        diff = np.abs(alnnum - alnnum[i])
        hit = np.sum((diff == 0)*1, axis=1)
        norm = np.sum((diff<25)*1, axis=1)
        seqid[i] = [hit[j]/norm[j] if norm[j] else 0 for j in range(aln.shape[0])]
    timelen = time.time() - ts
#    print(timelen)
    return seqid


aln_fn = sys.argv[1]
list_fn = sys.argv[2]
num = int(sys.argv[3])
numtot = int(sys.argv[4])
mutations_str = sys.argv[5]

names = []
with open(list_fn) as lf:
    for line in lf:
        names.append(line.split('|')[3])


wuhan_d = "2019-12-30"
wuhan_dt = datetime.datetime.strptime(wuhan_d, '%Y-%m-%d')
bin_size = 90
date_bins = []
this_dt = wuhan_dt
download_dt = datetime.datetime.strptime("2021-06-14", '%Y-%m-%d')
nseq_in_bins = {}
i = 0
color_list = ["rgba(255, 0, 0, 0.6)", "rgba(255, 165, 0, 0.6)", "rgba(255, 255, 0, 0.6)", "rgba(173, 255, 47, 0.6)", "rgba(0, 255, 0, 0.6)", "rgba(0, 255, 175, 0.6)", "rgba(0, 175, 255, 0.6)"]
colors = {}
legend = {}
while this_dt < download_dt:
    this_bin = (this_dt, this_dt + datetime.timedelta(days=bin_size))
    date_bins.append(this_bin)
    nseq_in_bins[this_bin] = 0 
    colors[this_bin] = color_list[i]
    legend[this_bin] = str(this_bin[0].date()) + " - " + str(this_bin[1].date())
    this_dt += datetime.timedelta(days=bin_size)
    i += 1

#print(date_bins)

lf = open("legend.txt", 'w')
lf.write("TREE_COLORS\nSEPARATOR TAB\nDATA\n")

nmax_in_bin = 100
nmaxtot = (nmax_in_bin+1)*len(nseq_in_bins)
nseqs = 0
tmpaln_fn = "tmpaln2.afa"
taf = open(tmpaln_fn, 'w')
seqs = set()
with open(aln_fn) as af:
    go = False
    iline = -1
    for line in af:
        iline += 1
#        if iline%10000 == 0:
#            print(nseq_in_bins)
        if nmaxtot == nseqs and not go:
            break
        if line.startswith(">"):
            go = False
            d = line.split('|')[2]
            if d[-2:] == '00':
                d = d[:-1] + '1'
            if d[-5:-3] == '00':
                d = d[:-4] + '1' + d[-3:]
            #print(d)
            dt = datetime.datetime.strptime(d, '%Y-%m-%d')
            for b in nseq_in_bins:
                if (dt >= b[0] and dt < b[1] and nseq_in_bins[b] <= nmax_in_bin and np.random.random() < 0.1) or iline == 0:
                    nseq_in_bins[b] += 1
                    nseqs += 1
                    go = True
                    if line.split('|')[3] not in names: ### TOGLI IL NOT
                        token = " #THIS#"
                    else:
                        token = ""
                    name = line.strip() + token + '\n'
                    curr_bin = b
                    break
        elif go and line.strip() not in seqs:
            taf.write(name)
            taf.write(line)
            seqs.add(line.strip())
            if "#THIS#" in name:
                token = "_THIS"
            else:
                token = ""
            lab = name.split('|')[3]+token
            lf.write("{0}\trange\t{1}\t{2}\n".format(lab, colors[b], legend[b]))
taf.close()
lf.close()

"""
aln, labels = parser_fasta(aln_fn, out_format="dict")
seqset = set()
newlabels = []
for name in labels:
    if aln[name] not in seqset:
        seqset.add(aln[name])
        newlabels.append(name)
print(len(labels), len(newlabels))
exit(1)
"""


"""
tmpaln_fn = "tmpaln.afa"
taf = open(tmpaln_fn, 'w')

seqset = set()
with open(aln_fn) as af:
    lastname = ""
    go = False
    for line in af:
        if line.startswith(">"):
            if np.random.random() < num/numtot:
                go = True
            else:
                go = False
            lastname = line.strip()
        elif go:
            if lastname.split('|')[3] not in names: ### TOGLI IL NOT
                token = " #THIS#"
            else:
                token = ""
            if line.strip() not in seqset: # NO DUPLICATE (ALIGNED) SEQUENCES
                #print(lastname.strip()+token)
                #print(line.strip())
                taf.write(lastname.strip()+token+"\n")
                taf.write(line.strip()+"\n")
                seqset.add(line.strip())
taf.close()
"""

aln, labels = parser_fasta(tmpaln_fn)
times_mx = np.zeros((len(labels), len(labels)))
times_mx2 = np.zeros((len(labels), len(labels)))
dt_arr, dt_arr2 = [], []
for ilabel, label in enumerate(labels):
    d = label.split('|')[2]
    if d[-2:] == '00':
        d = d[:-1] + '1'
    if d[-5:-3] == '00':
        d = d[:-4] + '1' + d[-3:]
    dt = datetime.datetime.strptime(d, '%Y-%m-%d')
    dt_arr.append(dt)
    dt_arr2.append(dt - wuhan_dt)
mindt = min(dt_arr)
maxdt = max(dt_arr)
dtspan = maxdt - mindt
dtspan = dtspan.days
dt_arr = np.array(dt_arr)
for i in range(len(labels)):
    for j in range(len(labels)):
        ddt = dt_arr[i] - dt_arr[j]
        ddt2 = dt_arr2[i].days + dt_arr2[j].days
        times_mx[i,j] = np.abs(ddt.days)/dtspan
        times_mx2[i,j] = np.abs(ddt2)

#print(times_mx)
"""
timealn, timelabels = {}, []
for name in labels:
    d = name.split('|')[2]
    if d[-1] == "0": # Unspecified day "00" is assumed to be first day of the month ("01")
        d = d[:-1] + "1"
    dt = datetime.datetime.strptime(d, '%Y-%m-%d')
    timelabels.append(name)
    timealn[name] = dt.date()
mintime = min([dt for dt in timealn.values()])
maxtime = max([dt for dt in timealn.values()])
timespan = maxtime-mintime
print(timespan.days)
# Divide in nbins bins
nbins = 5
binsinit_list = [mintime+datetime.timedelta(days=i) for i in range(0, timespan.days, int(timespan.days/(nbins+1)))][:-1]
listoflist = []
for ib in range(nbins):
    listoflist.append([])
alreadyin = set()

for ib in reversed(list(range(nbins))):
    for name in labels:
        if timealn[name] >= binsinit_list[ib] and name not in alreadyin:
            listoflist[ib].append(name)
            alreadyin.add(name)
minbinpop = min([len(l) for l in listoflist])
if minbinpop < 2:
    print("ERROR: every bin must have at least 2 sequences")
    exit(1)

print(minbinpop)
if minbinpop < int(num/nbins):
    print("WARNING: number of seqs reduced to {0} per bin (instead of {1})".format(minbinpop, int(num/nbins)))
    binpop = minbinpop
elif int(num/nbins) < 2:
    print("WARNING: number of seqs augmented to {0} per bin (instead of {1})".format(2, int(num/nbins)))
    binpop = 2
else:
    binpop = int(num/nbins)

rclist = []
for ib in reversed(list(range(nbins))):
    rclist += list(np.random.choice(listoflist[ib], size=binpop, replace=False))
#print(listoflist)
#print([len(l) for l in listoflist])
print(rclist)
"""

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
maxseqid = np.max(seqid)
seqid /= maxseqid
print(maxseqid)
print(seqid)
#dm = DistanceMatrix((1-seqid)*10, [x.split('|')[3]+"#"+x.strip()[-5:-1] for x in newlabels])
#dm = DistanceMatrix(((1-seqid)*50 + times_mx), [x.split('|')[3]+"#"+x.strip()[-5:-1] for x in newlabels])
print(times_mx)
plt.scatter(times_mx.flatten(), seqid.flatten())
plt.savefig("time_seqid.png")
plt.clf()
plt.scatter(times_mx2.flatten(), seqid.flatten())
plt.savefig("time2_seqid.png")
dm = DistanceMatrix(0.1*times_mx+seqid, newlabels)
nwk_t = tree.nj(dm, result_constructor=str)
print(Tree(nwk_t))
# We have to root the tree
t = Tree(nwk_t)
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
"""
nonames = []
rightmut = [b for a,j,b in mutations]
print(t)
t2 = t.copy('deepcopy')
for node in t.traverse():
    print(node.name, node.sankoff, node.sankoff == rightmut)
    if node.name in nonames:
        continue
    if node.sankoff == rightmut and len([x for x in node.iter_leaves()]) > 2:
        print(node.get_ascii(show_internal=True))
        nonames += [n.name for n in node.traverse()]
        t2.search_nodes(name=node.name)[0].detach()
print(t2)
"""

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
