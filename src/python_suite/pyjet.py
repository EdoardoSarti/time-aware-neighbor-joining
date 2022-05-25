from skbio import DistanceMatrix
from skbio import tree
import argparse
import matplotlib.pyplot as plt
from support import *
import sankoff

def create_sample(aln, names, seqids):
    M, L = aln.shape
    bins = [(0.2,0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 0.99)]
    name_bins = [[], [], [], []]
    for iname, name in enumerate(names):
        # Don't sample first row
        #if iname == 0:
        #    continue
        for ib, b in enumerate(bins):
            if seqids[iname] >= b[0] and seqids[iname] < b[1]:
                name_bins[ib].append(name)
                break

    sample_names, sample_rows = [], []
    for ib, b in enumerate(bins):
        errprint("[create_sample] REMARK: bin {0} has {1} sequences --> {2} sequences in sample tree".format(b, len(name_bins[ib]), int(len(name_bins[ib])**0.5)))
        sample_names += list(np.random.choice(name_bins[ib], int(len(name_bins[ib])**0.5), replace=False))
        sample_rows += [names.index(x) for x in sample_names]

    #sample_rows = [0] + sample_rows # First row is always included
    #sample_names = [sample_names[0]] + sample_names

    del_rows = [x for x in range(M) if x not in sample_rows]
    sample_aln = np.delete(aln, del_rows, 0)

    return sample_aln, sample_names


def jet_parse():
    parser = argparse.ArgumentParser()
    parser.add_argument("-A", "--alignment", type=str, help="input alignment file (fasta or stockholm formats)")
    args = parser.parse_args()
    return args


def jet():
    args = jet_parse()
    input_aln, labels = import_aln(args.alignment)  # imports fasta alignment (no identical seqnames allowed)
    refseq_aln, refseq_labels = adapt_aln_to_ref(input_aln, labels, ref=0)  # deletes columns where ref has gaps, puts ref in first place
    seqid_v = first_vs_all_seqid(refseq_aln)  # quick seqid between first ungapped seq and all the seqs (included itself)

    # Simplify labels, but keep them unique
    sp_refseq_labels = [x.split()[0].split('|')[0] for x in refseq_labels]
    # ...but still make them unique
    chd = {}
    for ilabel in range(len(sp_refseq_labels)):
            if sp_refseq_labels[ilabel] in chd:
                chd[sp_refseq_labels[ilabel]] += 1
                sp_refseq_labels[ilabel] += "_" + str(chd[sp_refseq_labels[ilabel]])
            else:
                chd[sp_refseq_labels[ilabel]] = 1

    alphabet = ['L','V','I','M','C','A','G','S','T','P','F','Y','W','E','D','N','Q','K','R','H','B','Z','X','-','.']
    sankoff_score_mx = sankoff.sankoff_subs_mx(alphabet) 
    init_d = sankoff.initialize_dictionary(alphabet)

    ctot = None
    nsamples = int(input_aln.shape[0]**0.5)
    for icount in range(nsamples):
        print("Create sample")
        sample_aln, sample_labels = create_sample(refseq_aln, refseq_labels, seqid_v)
        # list of simple labels
        sp_sample_labels = []
        for l in sample_labels:
            sp_sample_labels.append(sp_refseq_labels[refseq_labels.index(l)])

#        print(sample_aln, sample_labels)
        print("Calculate SeqID")
        seqid_mtx = []
        for i in range(len(sample_labels)):
            seqid_vi = first_vs_all_seqid(adapt_aln_to_ref(sample_aln, sample_labels, ref=i)[0])
            seqid_vi[[0, i]] = seqid_vi[[i, 0]]
            seqid_mtx.append(seqid_vi.copy())
        dist_mtx = np.ones(np.array(seqid_mtx).shape) - np.array(seqid_mtx)

        print("Create tree")
        dm = DistanceMatrix(dist_mtx, sp_sample_labels)
        nwk_t = tree.nj(dm, result_constructor=str)

        # We have to root the tree
        t = Tree(nwk_t)
        t.set_outgroup(t.get_midpoint_outgroup())
        nunks = 0
        for node in t.traverse():
            if not node.name.strip():
                node.name = 'UNK' + str(nunks).zfill(8)
                nunks += 1

        all_c = []

        print("Calculate Sankoff")
        leaves = {sp_sample_labels[i] : [(x.decode("utf-8"), init_d[x.decode("utf-8")]) for x in sample_aln[i]] for i in range(len(sample_labels))}
        t, ns, c = sankoff.sankoff(t, sankoff_score_mx, leaves, alphabet, jet=True)
        all_c.append(np.array(c))

        print("Calculate trace")
        c = np.array(sankoff.trace(t, alphabet, leaves))
        all_c.append(c)

        print("LENGTHS", input_aln.shape[1], refseq_aln.shape[1], sample_aln.shape[1], len(all_c[0]), len(all_c[1]))

        if type(ctot) == type(None):
            ctot = all_c
        else:
            ctot = [ctot[i] + all_c[i] for i in range(len(ctot))]

    ctot = [ctot[i]/nsamples for i in range(len(ctot))]
    ticklab = list(refseq_aln[0])
    ticks = [i for i in range(len(ticklab))]
    plt.xticks(ticks, ticklab)
    for c in ctot:
        plt.plot(c)
    plt.savefig("prova.png")


def fast_allvsall_seqid(aln, labels, alphabet, delchars={'-', '.', 'B', 'X', 'Z'}):
    real_alphabet = sorted(list(set(alphabet) - set(delchars)))
    conv = {c : real_alphabet.index(c) for c in real_alphabet}
    delscore = 50
    delstep = 25
    alnnum = np.zeros(aln.shape)
    for i in alnnum.shape[0]:
        for j in alnnum.shape[1]:
            if aln[i,j].decode('utf-8') in delchars:
                alnnum[i,j] = delscore
                delscore += delstep
            else:
                alnnum[i,j] = conv[aln[i,j].decode('utf-8')]
    print(alnnum)


def multigemme(aln_fn, mutations_str):
    # Parse mutations
    mutations = [(x[0], int(x[1:-1]), x[-1]) for x in sorted(mutations_str.split(','), key=lambda x: int(x[1:-1]))]  # Will fail if mutations_str not in the correct format <1-letter-aa><position><1-letter-aa>

    # Build aln mx
    input_aln, labels = import_aln(args.alignment)  # imports fasta alignment (no identical seqnames allowed)
    refseq_aln, refseq_labels = adapt_aln_to_ref(input_aln, labels, ref=0)  # deletes columns where ref has gaps, puts ref in first place
    seqid_v = first_vs_all_seqid(refseq_aln)  # quick seqid between first ungapped seq and all the seqs (included itself)

    # Simplify labels, but keep them unique
    sp_refseq_labels = [x.split()[0].split('|')[0] for x in refseq_labels]
    # ...but still make them unique
    chd = {}
    for ilabel in range(len(sp_refseq_labels)):
            if sp_refseq_labels[ilabel] in chd:
                chd[sp_refseq_labels[ilabel]] += 1
                sp_refseq_labels[ilabel] += "_" + str(chd[sp_refseq_labels[ilabel]])
            else:
                chd[sp_refseq_labels[ilabel]] = 1

    alphabet = ['L','V','I','M','C','A','G','S','T','P','F','Y','W','E','D','N','Q','K','R','H','B','Z','X','-','.']
    sankoff_score_mx = sankoff.sankoff_subs_mx(alphabet) 
    init_d = sankoff.initialize_dictionary(alphabet)

    # Sankoff
    print("Calculate SeqID")
    seqid_mtx = []
    for i in range(len(refseq_labels)):
        seqid_vi = first_vs_all_seqid(adapt_aln_to_ref(refseq_aln, sp_refseq_labels, ref=i)[0])
        seqid_vi[[0, i]] = seqid_vi[[i, 0]]
        seqid_mtx.append(seqid_vi.copy())
    dist_mtx = np.ones(np.array(seqid_mtx).shape) - np.array(seqid_mtx)

    print("Create tree")
    dm = DistanceMatrix(dist_mtx, sp_refseq_labels)
    nwk_t = tree.nj(dm, result_constructor=str)

    # We have to root the tree
    t = Tree(nwk_t)
    t.set_outgroup(t.get_midpoint_outgroup())
    nunks = 0
    for node in t.traverse():
        if not node.name.strip():
            node.name = 'UNK' + str(nunks).zfill(8)
            nunks += 1

    print("Calculate Sankoff")
    leaves = {sp_refseq_labels[i] : [(x.decode("utf-8"), init_d[x.decode("utf-8")]) for x in refseq_aln[i]] for i in range(len(refseq_labels))}
    t, ns = sankoff.sankoff(t, sankoff_score_mx, leaves, alphabet)

    # Take all the sequences in the aln that mutate ALL the specified positions
    t2 = t.copy('deepcopy')
    refleaf = t2.search_nodes(name=sp_refseq_labels[0])[0]
    t2.set_outgroup(refleaf)
    all_different = False
    reulting_trees = []
    for node in t2.traverse():
        cont = False
        for rt in reulting_trees:
            if node.name in [x.name for x in rt.traverse()]:
                cont = True
                break
        if cont:
            continue
        all_different = True
        for a, i, b in mutations:
            if node.sankoff[i-1] == a:
#                print("Node", node.name, a, i, b)
                all_different = False
                break
#        print(mutations)
#        print([node.sankoff[i-1] for a, i, b in mutations])
        if all_different:
            print("BECOMING ALL DIFFERENT AT", node.name)
            print(mutations)
#            print([node.sankoff[i-1] for a, i, b in mutations])
            print("BEFORE")
            print(t.get_ascii(show_internal=True))
            n = t.search_nodes(name=node.name)[0]
            found = False
            for c in n.children:
                if c.search_nodes(name=refleaf.name):
                    n.remove_child(c)
                    print("DELETE ITSELF", n.name)
                    newroot_name = None
                    if n != t:
                        newroot_name = n.up.name
                        print("NEW ROOT", newroot_name)
                    n.delete()
                    if type(newroot_name) != type(None):
                        t3 = t.copy('deepcopy')
                        newroot = t3.search_nodes(name=newroot_name)[0]
                        t3.set_outgroup(newroot)
                        reulting_trees.append(t3)
                    else:
                        reulting_trees.append(t.copy('deepcopy'))
                    found = True
                    break
            if not found: # it means it's coming from the root
                reulting_trees.append(n.copy('deepcopy'))
            else:
                print(reulting_trees[-1].get_ascii(show_internal=True))
                break
            print(reulting_trees[-1].get_ascii(show_internal=True))

    # Stat: analysis of subtrees and how many right mutations they conserve (i.e. internal node has right mutation, and there are at least 3 leaves that are connected to that node by always conserving the residue)
    # how many mutations vs. how big, how deep
    right_muts = {}
    for n_of_right_muts in range(1,len(mutations)+1):
        print(n_of_right_muts, "CORRECT MUTATIONS")
        right_muts_trees = []
        for rt in reulting_trees:
            for node in rt.traverse():
                cont = False
                for rmt in right_muts_trees:
                    if node.name in [x.name for x in rmt.traverse()]:
                        cont = True
                        break
                if cont:
                    continue
                print(node.name)
                print(mutations)
                print([node.sankoff[i-1] for a, i, b in mutations])
                if node.name not in right_muts:
                    right_muts[node.name] = []
                    for a, i, b in mutations:
                        if node.sankoff[i-1] == b:
                            right_muts[node.name].append((a,i,b))
                if len(right_muts[node.name]) == n_of_right_muts:
                    right_muts_trees.append(node.copy('deepcopy'))
                    print(node.get_ascii(show_internal=True))


    # The chosen tree is the one that maximises the graph in the analysis

    # On that aln, make a new NJ tree and Sankoff only of the mutation positions

    # Compare the two trees and calculate how many changes there are throughout those trees in those positions

if __name__ == "__main__" :
#    jet()
    parser = argparse.ArgumentParser()
    parser.add_argument("-A", "--alignment", type=str, help="input alignment file (fasta or stockholm formats, first sequence is reference, must be without gaps)")
    parser.add_argument("-m", "--mutations", type=str, help="mutations csv string (X#Y format, numbering starts from 1, X must reflect residue in first sequence of aln)")
    args = parser.parse_args()
    multigemme(args.alignment, args.mutations)
