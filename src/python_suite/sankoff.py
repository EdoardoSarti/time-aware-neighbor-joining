# -*- coding: utf-8 -*-
"""
Created on Jun 01 2021

@author: sarti
"""
from support import *


def sankoff_subs_mx(alphabet, equivalence_class=['B','Z','X','.','-']):
    mx_levels = [[['L'],['V'],['I'],['M'],['C'],['A'],['G'],['S'],['T'],['P'],['F'],['Y'],['W'],['E'],['D'],['N'],['Q'],['K'],['R'],['H']],
        [['L','V','I','M'],['C'],['A'],['G'],['S'],['T'],['P'],['F','Y'],['W'],['E'],['D'],['N'],['Q'],['K','R'],['H']],
        [['L','V','I','M'],['C'],['A'],['G'],['S','T'],['P'],['F','Y','W'],['E','D','N','Q'],['K','R'],['H']],
        [['L','V','I','M','C'],['A','G','S','T','P'],['F','Y','W'],['E','D','N','Q','K','R','H']],
        [['L','V','I','M','C','A','G','S','T','P','F','Y','W'],['E','D','N','Q','K','R','H']],
        [['L','V','I','M','C','A','G','S','T','P','F','Y','W','E','D','N','Q','K','R','H']]]
    mx_d = {}
    for a in sorted(mx_levels[-1][0]):
        for b in sorted(mx_levels[-1][0]):
            brk = False
            for ilev, lev in enumerate(mx_levels):
                for l in lev:
                    if a in l and b in l:
                        mx_d[(a,b)] = ilev*2
                        brk = True
                        break
                if brk:
                    break

    L = len(alphabet)
    subscore = np.ones((L,L))*INF
    for ia, a in enumerate(alphabet):
        for ib, b in enumerate(alphabet):
            if (a,b) in mx_d:
                subscore[ia,ib] = mx_d[(a,b)]
            elif a == b:
                subscore[ia,ib] = 0
            elif (a in equivalence_class and b not in equivalence_class):
                subscore[ia, ib] = 1

    return subscore


def initialize_dictionary(alphabet, ones=False):
    # Dictionary to store node parcimony scores, initialized by leaves scores
    l = len(alphabet)
    if ones:
        d = {c : np.array([0 if k != i else 1 for k in range(l)])
            for i, c in enumerate(alphabet)}
    else:
        d = {c : np.array([float(INF) if k != i else 0 for k in range(l)]) 
            for i, c in enumerate(alphabet)}
    return d


def pjoin(node_1, node_2, score_mx, idx):
    parent = node_1.up
#    print("THIS PARENT", parent, parent.sankoff_score)
    L = score_mx.shape[0]

    if parent != node_2.up:
        print("ERROR: the two nodes are not siblings", node_1.name, node_2.name)
        exit(1)

    allminsc = INF
    for i in range(L): # on the parent
        minsc = 0
#        parent.traces.append([])
        traces = {}
        for thisnode in [node_1, node_2]:
            sumv = thisnode.sankoff_score[idx] + score_mx[:,i]
#            print("parent_char:", i, thisnode.name, sumv)
            argminsc = np.argmin(sumv)
            minsc += sumv[argminsc]
            traces[thisnode.name] = argminsc
        parent.sankoff_score[idx][i] = minsc
        if minsc < allminsc:
            allminsc = minsc
            parent.traces = traces
#        print("PARENT UPDATE", parent.sankoff_score[idx], parent.traces)
         

def traverse_order(tree, only_topology=False):
    if not only_topology:
        # Trees maintain distances, and levels take them into consideation. Typically there is 1 node each level
        # WARNING: leaves are not considered
        distances = [(node, tree.get_distance(node)) for node in tree.traverse() if not node.is_leaf()]
        distances = sorted(distances, key=lambda x : x[1])
#        print("DISTANCES", distances)

        prev_d = -9999
        traverse = []
        level = []
        for x in distances:
            level.append(x[0])
            if x[1] != prev_d:
                traverse.append(level)
                level = []
    else:
        # Only topology is preserved
        # WARNING: leaves are not considered
        stop = False
        this_level = [tree]
        next_level = []
        traverse = []
        traverse.append(this_level)
        while True:
            for node in this_level:
                next_level += [x for x in node.children if not x.is_leaf()] 
            if next_level:
                this_level = next_level
                traverse.append(this_level)
#                print(len(this_level), this_level)
                next_level = []
            else:
                break

    # Add names
    traverse_names = []
    for l in traverse:
        ll = []
        for n in l:
            ll.append(n.name)
        traverse_names.append(ll)

    return traverse, traverse_names


def trace(tree, alphabet, leaves_d, discard_char=['-','X','.', 'Z', 'B']):
#    print("LEAVES", [leaves_d[n.name][0][0] for n in tree.iter_leaves()])

    # Order
    order, order_names = traverse_order(tree)

    # Do not consider these characters in the alphabet
    reducedalph = sorted(list(set(alphabet) - set(discard_char)))
    init_d = initialize_dictionary(reducedalph, ones=True)

    Ll = len(leaves_d[[x for x in leaves_d][0]])
    totcons, totcons_d = [], {i : 0 for i in range(Ll)}

    # Build an order (from leaves to root, level by level, pairs of children)
    order, order_names = traverse_order(tree)

    # Initialize leaves
    for leaf in tree.iter_leaves():
        leaf.add_features(c=np.array([init_d[leaves_d[leaf.name][x][0]] if leaves_d[leaf.name][x][0] in reducedalph else np.zeros(len(reducedalph)) for x in range(Ll)]))
    for ilevel, level in [(x,y) for x,y in enumerate(order)][::-1]:
        for node in level:
            m1, m2 = [chi.c for chi in node.children]
            node.add_features(c=m1+m2)
                
    respos = {i for i in range(Ll)}
    respos_l = np.zeros(Ll)
    for ilevel, level in [(x,y) for x,y in enumerate(order)]:
        for node in level:
            # Subtrees with 2 leaves or less are not considered for conservation
            liter = list(node.iter_leaves())
            lenliter = len(liter)
            if lenliter < 3:
                continue
            for i in list(respos):
                if np.max(node.c[i]/lenliter) > 0.9:
                    respos_l[i] += 1
                if respos_l[i] == 2:
                    totcons_d[i] = (len(order)-ilevel)/len(order)
                    respos.discard(i)
    totcons = [totcons_d[i] for i in range(Ll)]
#    print(totcons)
    return totcons


def read_score(score_fn):
    # Import score file: 3 columns (row_index, column_index, value)
    score_d = {}
    alphabet = []
    with open(score_fn) as sf:
        for line in sf:
            if (not line.strip()) or line.strip().startswith("#"):
                continue
            fields = line.split()
            if fields[0] not in score_d:
                score_d[fields[0]] = {}
                alphabet.append(fields[0])
            score_d[fields[0]][fields[1]] = float(fields[2])
    alphabet = sorted(alphabet)
    A = len(alphabet)
    score = np.zeros((A,A))
    for i1, c1 in enumerate(alphabet):
        for i2, c2 in enumerate(alphabet):
            score[i1,i2] = score_d[c1][c2] # numpy square score matrix ordered as alphabet
    return score, alphabet


def build_sankoff(tree_fn, score_fn, leaves_fn):
    # Import ETE3 tree (SUPPOSES ALL NODES HAVE A NAME, except root)
    tree = Tree(tree_fn)
    N = len([x for x in tree.traverse()])
    names_sank = {}
    if not tree.name.strip():
        tree.name = "ROOT"

    score, alphabet = read_score(score_fn)

    # Import leaf information file: 2 columns (leaf_name, character)
    leaves_d = {}
    Ll = 1
    init_d = initialize_dictionary(alphabet)  # defines leaf vectors for each character of the alphabet
    with open(leaves_fn) as lf:
        for line in lf:
            if (not line.strip()) or line.strip().startswith("#"):
                continue
            fields = line.split()
            l = list(fields[1])
            Ll = len(l)
            leaves_d[fields[0]] = [(x, init_d[x]) for x in l] # dict leaf_name : (chacacter, leaf_vector(character))

    # Check that leaves_fn had all and only the leaf names

    return sankoff(tree, score, leaves_d, alphabet)


def sankoff(tree, score, leaves_d, alphabet, jet=False, jet_discard=['-','X','.']):
#    print("LEAVES", [leaves_d[n.name][0][0] for n in tree.iter_leaves()])
    A = len(alphabet)
    INF_V = np.ones(A)*INF
    Ll = len(leaves_d[[x for x in leaves_d][0]])
    N = len([x for x in tree.traverse()])
    names_sank = {node.name : [] for node in tree.traverse()}
    if jet:
        totcons = []

    # Build an order (from leaves to root, level by level, pairs of children)
    order, order_names = traverse_order(tree)

    # Initialize Sankoff features
    for node in tree.traverse():
        if node.is_leaf():
            node.add_features(sankoff=[leaves_d[node.name][i][0] for i in range(Ll)], sankoff_score=np.array([leaves_d[node.name][i][1] for i in range(Ll)]), traces={})
        else:
            node.add_features(sankoff=['0']*Ll, sankoff_score=[np.copy(INF_V)]*Ll, traces={})


    # Loop on the position in the leaves' sequences
    for il in range(Ll):
        for ilevel, level in [(x,y) for x,y in enumerate(list(reversed(order)))]:
            for inode, node in enumerate(level):
                node.traces = {}
#        print("LEAVES AGAIN", [n.sankoff[il] for n in tree.iter_leaves()])

        # Traverses the tree, level by level, starting from the leaves and going up to the root-1
        for ilevel, level in [(x,y) for x,y in enumerate(list(reversed(order)))]:
            for inode, node in enumerate(level):
                #print(node.name, [x.name for x in node.children])
                node1, node2 = node.children
                pjoin(node1, node2, score, il)

        # Traceback    
        order[0][0].sankoff[il] = alphabet[np.argmin(order[0][0].sankoff_score[il])]
        names_sank[order[0][0].name].append(order[0][0].sankoff[il])
#        print("ROOT", order[0][0].sankoff)
        for ilevel, level in [(x,y) for x,y in enumerate(order)][1:]:
#            print("LEVEL", ilevel)
            for inode, node in enumerate(level):
                node.sankoff[il] = alphabet[node.up.traces[node.name]]
                names_sank[node.name].append(node.sankoff[il])
#                print(node.sankoff[il], node.up.traces, node.sankoff_score[il])

        if jet:
            # Find % conserved subtrees
            consST = 0
            cons = 0
            brk = False
            jetalph = list(set(alphabet) - set(jet_discard))
            for ilevel, level in [(x,y) for x,y in enumerate(order)]:
                for node in level:
                    # Subtrees with 2 leaves or less are not considered for conservation
                    if len(list(node.iter_leaves())) < 3 or brk:
                        continue
                    stats = {a : 0 for a in alphabet}
#                    print([(x,y.name) for x,y in enumerate(node.iter_leaves())])
                    for ileaf, leaf in enumerate(node.iter_leaves()):
                        stats[leaf.sankoff[il]] += 1
                    norm = ileaf+1
                    for a in jetalph:
#                        print(node.name, a, stats[a]/norm)
                        if stats[a]/norm > 0.9:
#                            print("LEVEL FOUND:", ilevel, "LEAVES:", stats[a])
                            consST += 1
                            if consST == 2:
#                                print(tree)
                                cons = (len(order)-ilevel)/len(order)
#                                print("CONSERVATION:", cons)
                                brk = True
            totcons.append(cons)

    for ilevel, level in [(x,y) for x,y in enumerate(order)]:
        for inode, node in enumerate(level):
            tree.search_nodes(name=node.name)[0].add_features(sankoff=node.sankoff)

    if jet:
        return tree, names_sank, totcons

    return tree, names_sank


# =============================================================================
# Test
# =============================================================================
#


def test():
    tree_fn = sys.argv[1]
    tree = Tree(tree_fn)
    nunks = 0
    for node in tree.traverse():
        if not node.name.strip():
            node.name = 'UNK' + str(nunks).zfill(8)
            nunks += 1
    tree.write(features=['name', 'sankoff'], format=0, outfile='inttree.txt')
    
    score_fn = 'score_try.txt'#sys.argv[2]
    leaves_fn = 'leaves_try.txt'#sys.argv[3]
    
    # Alphabet:
    alph = ['A','C','G','T']
    # Score:
    with open(score_fn, 'w') as sf:
        for i in range(len(alph)):
            for j in range(len(alph)):
                sf.write("{0}\t{1}\t{2}\n".format(alph[i], alph[j], 0 if i==j else np.random.uniform(low=0.5, high=2)))
    # Leaves:
    nl = []
#    print(tree)
    with open(leaves_fn, 'w') as lf:
        for leaf in tree.iter_leaves():
            rc = []
            for i in range(350):
                rc.append(np.random.choice(alph))
            nl.append((leaf.name, rc))
            lf.write("{0}\t{1}\n".format(leaf.name, "".join(rc)))
    
    t, ns = build_sankoff('inttree.txt', score_fn, leaves_fn)
    t.write(features=['name', 'sankoff'], format=0, outfile="out_tree_try.nw")
    with open('outseqs', 'w') as of:
        for node in t.traverse():
            of.write(">"+node.name+"\n")
            of.write("".join(ns[node.name])+"\n")


if __name__ == "__main__":
    test()
