from skbio import DistanceMatrix
from skbio import tree
from ete3 import Tree
import numpy as np
import re

# PERCHE DISTANZE NEGATIVE???

def internal_nj(dist_mx, labels):
    # Run checks on distance matrix and label association
    n = len(labels)
    if dist_mx.shape != (n,n):
        print("[internal_nj] ERROR: distance MX shape {0}, number of labels {1}".format(dist_mx.shape, n))
        exit(1)

    #  Symmetry
    for i in range(n):
        for j in range(n):
            if dist_mx[i,j] != dist_mx[j,i]:
                print("DISTANCES NOT SYMMETRIC", i, j)

    #  Triangle inequality
    eps = 0.001
    for i in range(n):
        for j in range(i):
            if not all(dist_mx[i,j] <= dist_mx[i,:] + dist_mx[:,j] + eps):
                for k in range(n):
                    print(dist_mx[i,j] <= dist_mx[i,k] + dist_mx[k,j] + eps, dist_mx[i,j], dist_mx[i,k], dist_mx[k,j])
                exit(1)
            

    # Initialize result string
    res = ""

    # Create working distance matrix
    row_labels = labels[:] # Labels of the plain distance / NJ matrix
    col_labels = labels[:] # Labels of the workall matrix (leaves are kept)
    work_mx = dist_mx.copy()  # The distance matrix whose columns never get reduced 
    parent = {i : i for i in range(work_mx.shape[1])} # Parent of each column label of work_mx

    newid = -1
    PRECISION = 5
    THR = 0.1**PRECISION
    while work_mx.shape[0] > 2:	# Rows get reduced as in plain NJ
        n, m = work_mx.shape

        # Generalized Neighbor Joining matrix (NJ plus rectangular complement) 
        si = np.sum(work_mx[:,:n], axis=1) # Sum on columns, leave rows as variable
        gnj_mx = 99999.*np.ones(work_mx.shape)
        for i in range(work_mx.shape[0]):
            # Plain NJ matrix
            for j in range(work_mx.shape[0]):
                if i == j:
                    continue
                gnj_mx[i,j] = (n-2)*work_mx[i,j] - si[i] - si[j] 
            # Complement
            for j in range((work_mx.shape[0], work_mx.shape[1]):
                gnj_mx[i,j] = 2*(n-2)*work_mx[i,j] - 2*si[parent[j]] - 2*(n-2)*work_mx[i,parent[j]]

        # Find minimum of the current GNJ matrix
        min_ind = np.unravel_index(np.argmin(gnj_mx, axis=None), gnj_mx.shape)
        search_for_accessible_leaves = False
        if min_ind[1] < gnj_mx.shape[0]:
            # Plain NJ criterion
            new_move = False
            newd_1 = 0.5*work_mx[min_ind] + 1/(2*(n-2))*(si[min_ind[0]] - sj[min_ind[1]])
            newd_2 = work_mx[min_ind] - newd_1
            if not ((newd_1 <= THR and newd_2 > 2*THR) or (newd_2 <= THR and newd_1 > 2*THR)):
                search_for_accessible_leaves = True
        else:
            # Creation of a degeneracy
            new_move = True
            newd_1 = work_mx[min_ind]
            newd_2 = 0	# At distance 0 from the old leaf (now root of a new subtree)
        newd = np.array([newd_1, newd_2])


        # Check if new node is degenerate and write new part of text to add to tree
        degenerate = True
        if new_move:
            tree_txt = "({0}:{1:.5f},{2}:{3:.5f})".format(row_labels[min_ind[0]], newd_1, col_labels[min_ind[1]], newd_2)
            degenerate = False
        elif row_labels[min_ind[0]].split(")")[-1][:3] == "TMP" and col_labels[min_ind[1]].split(")")[-1][:3] == "TMP":
            t1 = Tree(row_labels[min_ind[0]]+";", format=1)
            t2 = Tree(col_labels[min_ind[1]]+";", format=1)
            new_child1 = "{0}:{1:.5f}".format(t1.name, newd_1)
            new_child1 = "{0}:{1:.5f}".format(t2.name, newd_2)
            if t1.get_distance(t1.children[0]) > t1.get_distance(t1.children[1]):
                new_grandchild1 = "{0}".format(t1.children[0].write(format=1)[:-1])
            else:
                new_grandchild1 = "{0}".format(t1.children[1].write(format=1)[:-1])
            if t2.get_distance(t2.children[0]) > t2.get_distance(t2.children[1]):
                new_grandchild2 = "{0}".format(t2.children[0].write(format=1)[:-1])
            else:
                new_grandchild2 = "{0}".format(t2.children[1].write(format=1)[:-1])
            tree_txt = "(({0}){1}, ({2}){3})".format(new_grandchild1, new_child1, new_grandchild2, new_child2)  # Join the two nodes 
            degenerate = False
        elif row_labels[min_ind[0]].split(")")[-1][:3] == "TMP":
            t1 = Tree(row_labels[min_ind[0]]+";", format=1)
            for c in t1.children:
                if t1.get_distance(c) <= THR:
                    new_parent = c.name
                else:
                    new_child1 = "{0}".format(c.write(format=1)[:-1])
                    new_child2 = "{0}:{1:.5f}".format(col_labels[min_ind[1]], newd_1+newd_2)
            wdim, ldim = 0, 1
            tree_txt = "({0},{1}){2}".format(new_child1, new_child2, new_parent)  # Join the two nodes
        elif col_labels[min_ind[1]].split(")")[-1][:3] == "TMP":
            t1 = Tree(col_labels[min_ind[1]]+";", format=1)
            for c in t1.children:
                if t1.get_distance(c) <= THR:
                    new_parent = c.name
                else:
                    new_child1 = "{0}:{1:.5f}".format(row_labels[min_ind[0]], newd_1+newd_2)
                    new_child2 = "{0}".format(c.write(format=1)[:-1])
            wdim, ldim = 1, 0
            tree_txt = "({0},{1}){2}".format(new_child1, new_child2, new_parent)  # Join the two nodes
        else:
            degenerate = False
            new_parent = ""
            new_child1 = "{0}:{1:.5f}".format(row_labels[min_ind[0]], max(THR, newd_1))
            new_child2 = "{0}:{1:.5f}".format(col_labels[min_ind[1]], max(THR, newd_2))
            tree_txt = "({0},{1}){2}".format(new_child1, new_child2, new_parent)  # Join the two nodes

        

#        print("TREETEXT", tree_txt)
        print(Tree(tree_txt+";", format=1).get_ascii(show_internal=True))
        print("NEW DISTANCES", newd)



        # Search a possible leaf to associate to the internal node
        assigned = None
        if search_for_accessible_leaves and not new_move:
            accessible_leaves = [x for x in labels if ("("+x+":" not in tree_txt) and (","+x+":" not in tree_txt) and (")"+x+":" not in tree_txt) and (x in col_labels)]  # internal nodes can be filled only by a subset of leaves
    
            # Is there any ORIGINAL LEAF NOT IN THE TWO SUBTREES that respects these distances? Check
            # Calculate all distances between ORIGINAL LEAVES MINUS LEAVES IN SUBTREES and the new node
            # Use this array in the next loop, instead of work_mx
            t = []
            for label in accessible_leaves:
                icol = col_labels.index(label)
                d = np.array([work_mx[min_ind[0], icol], workall_mx[min_ind[1], icol]]) # search_for_accessible_leaves tag guarantees min_ind[1] < work_mx.shape[0]
                t.append((np.linalg.norm(d - newd)/np.linalg.norm(newd), icol))  # Distances from the rest of the tree and difference between branch lengths
            t = sorted(t, key= lambda x: x[0])  # The difference between the branch length must be minimized. Optimal: 0, Threshold: 1 (difference must be less than original branch length)
     
            # Stub for more complex criterion of choice
            if t[0][0] < 1:
                assigned = t[0]


        # Manage matrices and label lists according to one of the five possibilities:
        if degenerate:
            # 1) Previously declared degeneracy to be resolved
            print("Resolution of one degeneracy")
            # This is the only case where also the second distance matrix loses a row/column
            # It's due to the fact that this move equilibrates the creation of a placeholder
            # (see first 2 options)
            wlabel, llabel = col_labels[min_ind[wdim]], col_labels[min_ind[ldim]]
            
            u_idxall_v = treeall_labels.index(tree_labels[vincente])
            u_idxall_p = treeall_labels.index(tree_labels[perdente])

            tree_labels[vincente] = tree_txt

            # Delete row/column/label in "plain" distance matrix

            # Rename the row/column in second matrix
            treeall_labels[u_idxall_v] = tree_txt

            nope = True
            for ipair, pair in enumerate(nojoin):
                if vincente == pair[1]:
                    nope = False
#                    print("REPLACE IN ", tree_labels[pair[0]], "THIS", treeall_labels[u_idxall_v].split(")")[-1], "WITH", tree_labels[vincente])
                    treeall_labels[treeall_labels.index(tree_labels[pair[0]])] = treeall_labels[treeall_labels.index(tree_labels[pair[0]])].replace(treeall_labels[u_idxall_v].split(")")[-1], tree_labels[vincente])
                    tree_labels[pair[0]] = tree_labels[pair[0]].replace(treeall_labels[u_idxall_v].split(")")[-1], tree_labels[vincente])
                    work_mx = np.delete(work_mx, [perdente, pair[1]], 0)
                    work_mx = np.delete(work_mx, [perdente, pair[1]], 1)
                    workall_mx = np.delete(workall_mx, [treeall_labels.index(tree_labels[pair[1]])], 0)
                    workall_mx = np.delete(workall_mx, [treeall_labels.index(tree_labels[pair[1]])], 1)
                    del treeall_labels[treeall_labels.index(tree_labels[pair[1]])]
                    newnojoin = []
                    for ppair in nojoin:
                        if ppair == pair:
                            continue
                        p1, p2 = ppair
                        if p1 > perdente:
                            p1 -= 1
                        if p2 > perdente:
                            p2 -= 1
                        if p1 > pair[1]:
                            p1 -= 1
                        if p2 > pair[1]:
                            p2 -= 1
                        newnojoin.append((p1, p2))
                    for i in sorted([perdente, pair[1]], reverse=True):
                        del tree_labels[i]
                    break
            nojoin = newnojoin
            if nope:
#                print("NOPE")
                del tree_labels[perdente]
                newnojoin = []
                for pair in nojoin:
                    p1, p2 = pair
                    if p1 > perdente:
                        p1 -= 1
                    if p2 > perdente:
                        p2 -= 1
                    newnojoin.append((p1, p2))
                nojoin = newnojoin
                work_mx = np.delete(work_mx, [perdente], 0)
                work_mx = np.delete(work_mx, [perdente], 1)
                

            # Delete row/column/label in second distance matrix
#            del treeall_labels[u_idxall_p]
#            workall_mx = np.delete(workall_mx, [u_idxall_p], 0)
#            workall_mx = np.delete(workall_mx, [u_idxall_p], 1) 
        elif (tree_labels[min_ind_all[0]] in labels and (newd_1 <= THR and newd_2 > 2*THR)):
            # 2) New degeneracy (left side)
            # The move is analogous to the classical NJ, but the new node IS the leaf with distance <= Threshold
            # - The row/column with distance > Threshold + unit is renamed in order for it to be the new node
            # - Deletes the two row/column in 'plain' distances, associated to the leaf with distance <= Threshold
            # - A corresponding renaming happens in the 'second' distance matrix
            # NOTA BENE: THE NEW NODE **IS** ONE OF THE OLD LEAVES -> NO NJ DISTANCE RECALCULATION IS MADE

            print("Creation left-degeneracy")
#            print("NEWD", treeall_labels[min_ind_all[0]], newd_1, treeall_labels[min_ind_all[1]], newd_2)
            u_idxall = treeall_labels.index(treeall_labels[min_ind_all[0]])
            # Rename the row/column (internal)
            newtxt = "TMP{0}".format(treeall_labels[min_ind_all[0]])
            treeall_labels[min_ind_all[0]] = tree_txt + newtxt

            # Delete row/column/label in "plain" distance matrix
            del treeall_labels[min_ind_all[1]]
            newnojoin = []
            for pair in nojoin:
                p1, p2 = pair
                if p1 > min_ind_all[1]:
                    p1 -= 1
                if p2 > min_ind_all[1]:
                    p2 -= 1
                newnojoin.append((p1, p2))
            nojoin = newnojoin
            work_mx = np.delete(work_mx, [min_ind_all[1]], 0)  # Delete the row
            work_mx = np.delete(work_mx, [min_ind_all[1]], 1)  # Delete the column

            # Rename the row/column in second matrix
            treeall_labels[u_idxall] = tree_txt + newtxt

#            print(tree_labels, treeall_labels)
        elif new_move or (tree_labels[min_ind[1]] in labels and (newd_1 > 2*THR and newd_2 <= THR)):
#        elif (tree_labels[min_ind[1]] in labels and (newd_1 > 2*THR and newd_2 <= THR)):
            # 3) New degeneracy (right side)
            print("Creation right-degeneracy")
            if new_move:
                # Cancella left vector, fai ricomparire right vector in work_mx e chiamalo (left:D, right:0)TMPright
                u_idx_1 = tree_labels.index(treeall_labels[min_ind_all[0]]) # Must be there
                newvec = np.array([workall_mx[min_ind_all[1],i] for i in range(workall_mx.shape[0]) if treeall_labels[i] in tree_labels])
                work_mx[u_idx_1,:] = newvec 
                work_mx[:,u_idx_1] = newvec
                tree_labels[u_idx_1] = tree_txt + "TMP{0}".format(treeall_labels[min_ind_all[1]])

                # Prevent work_mx to join the two severed branches (they will be joined again once the TMP has been resolved)
                for il, l in enumerate(tree_labels):
                    if treeall_labels[min_ind_all[1]] in l and treeall_labels[min_ind_all[1]] != l:
                        nojoin.append((il, u_idx_1))

                # In workall, right label becomes new label (vectors are maintained)
                treeall_labels[min_ind_all[1]] = tree_txt + "TMP{0}".format(treeall_labels[min_ind_all[1]])
            else:    
#                print("NEWD", tree_labels[min_ind[1]], newd_2, tree_labels[min_ind[0]], newd_1)
                u_idxall = treeall_labels.index(tree_labels[min_ind[1]])
                # Rename the row/column (internal)
                newtxt = "TMP{0}".format(tree_labels[min_ind[1]])
                tree_labels[min_ind[1]] = tree_txt + newtxt

                # Delete row/column/label in "plain" distance matrix
                del tree_labels[min_ind[0]]
                newnojoin = []
                for pair in nojoin:
                    p1, p2 = pair
                    if p1 > min_ind[0]:
                        p1 -= 1
                    if p2 > min_ind[0]:
                        p2 -= 1
                    newnojoin.append((p1, p2))
                nojoin = newnojoin
                work_mx = np.delete(work_mx, [min_ind[0]], 0)  # Delete the row
                work_mx = np.delete(work_mx, [min_ind[0]], 1)  # Delete the column

                # Rename the row/column in second matrix
                treeall_labels[u_idxall] = tree_txt + newtxt
        elif type(assigned) != type(None): 
            # 4) Leaf assigned to be the new internal node

            print("Internal node assignment")
            # No row/column to add: just one to rename and one to delete
            u_idxall = assigned[1]
#            print("LEAF FOUND", treeall_labels[u_idxall])
            if assigned[2] != -1:
                # Rename the row/column (external)
                tree_labels[assigned[2]] = tree_txt + treeall_labels[u_idxall]

                # Delete rows/columns/labels in "plain" distance matrix
                for i in sorted(min_ind, reverse=True):
                    del tree_labels[i]
                newnojoin = []
                for pair in nojoin:
                    p1, p2 = pair
                    if p1 > min_ind[0]:
                        p1 -= 1
                    if p2 > min_ind[0]:
                        p2 -= 1
                    if p1 > min_ind[1]:
                        p1 -= 1
                    if p2 > min_ind[1]:
                        p2 -= 1
                    newnojoin.append((p1, p2))
                nojoin = newnojoin
                work_mx = np.delete(work_mx, min_ind, 0)
                work_mx = np.delete(work_mx, min_ind, 1)

                # Rename the row/column in second matrix
                treeall_labels[u_idxall] = tree_txt + treeall_labels[u_idxall]
            else:
                # This means that that leaf has already been grouped, thus the column in the NJ matrix has gone
                # We thus have to retrieve the node that contains the information of the wanted leaf
#                print("LOOK FOR", treeall_labels[u_idxall])
#                print("IN", tree_labels)
                ii = [treeall_labels[u_idxall] in x for x in tree_labels].index(True) # Look for the only label in the plain NJ tree containing the wanted label
#                print("FOUND AT", ii)
#                print(tree_labels[ii])

                # Rename the row/column
#                print("BEF", tree_labels[ii])
                tree_labels_tmp = tree_labels[ii].replace(treeall_labels[u_idxall], tree_txt + treeall_labels[u_idxall])
                tree_labels[ii] = tree_labels_tmp
#                print("AFT", tree_labels[ii])

                # Delete rows/columns/labels in "plain" distance matrix
                for i in sorted(min_ind, reverse=True):
#                    print("DEL", i, tree_labels[i])
                    del tree_labels[i]
                newnojoin = []
                for pair in nojoin:
                    p1, p2 = pair
                    if p1 > min_ind[0]:
                        p1 -= 1
                    if p2 > min_ind[0]:
                        p2 -= 1
                    if p1 > min_ind[1]:
                        p1 -= 1
                    if p2 > min_ind[1]:
                        p2 -= 1
                    newnojoin.append((p1, p2))
                nojoin = newnojoin
                work_mx = np.delete(work_mx, min_ind, 0)
                work_mx = np.delete(work_mx, min_ind, 1)

                # Rename the row/column in second matrix
                treeall_labels[u_idxall] = tree_labels_tmp#treeall_labels[u_idxall].replace(treeall_labels[u_idxall], tree_txt + treeall_labels[u_idxall])

            for i in sorted(min_ind_all, reverse=True):
                if treeall_labels[i] not in labels:
                    del treeall_labels[i]
                    workall_mx = np.delete(workall_mx, [i], 0)
                    workall_mx = np.delete(workall_mx, [i], 1)
        else:
            # 5) Classical neighbor joining
            print("Creation of new internal node (classic NJ)")
            # Plain NJ option
            # - Adds one row/column in 'plain' distances, associated to a new node, with the classical NJ formula
            # - Deletes the two rows/columns in 'plain' distances, associated to the two joined subtrees
            # - Adds one corresponding row/column in 'second' distance matrix, with the classical NJ formula (yet calculated on the 'second' matrix instead)

            # Distances new row/column
            u_col = np.zeros(work_mx.shape[0])
            for i in range(work_mx.shape[0]):
                u_col[i] = 0.5*(work_mx[min_ind[0], i] + work_mx[i, min_ind[1]] - work_mx[min_ind[0], min_ind[1]]) # Classical NJ formula

            # Distances new row/column in second matrix
            u_row = np.zeros(work_mx.shape[1])
            u_row[:work_mx.shape[0]] = u_col
            for i in range(work_mx.shape[0], work_mx.shape[1]):
                # dist(new, oldleaf) = dist(new,parent(oldleaf)) + dist(parent(oldleaf), oldleaf) = NJdist(new,parent(oldleaf)) + dist(parent(oldleaf), oldleaf)
                u_row[i] = 0.5*(work_mx[min_ind[0], parent[i]] + work_mx[parent[i], min_ind[1]] - work_mx[min_ind[0], min_ind[1]]) + work_mx[parent[i], i]

            # ARRIVATO QUI
            # DEVO ANCORA CAMBIARE OPZIONI 1-4


            # Add row/column to matrices
            work_mx = np.vstack((work_mx, u_dist))
            workall_mx = np.vstack((workall_mx, uall_dist))
            u_dist = np.hstack((u_dist, [0]))
            work_mx = np.column_stack((work_mx, u_dist))
            uall_dist = np.hstack((uall_dist, [0]))
            workall_mx = np.column_stack((workall_mx, uall_dist))

            # Update labels of matrices
            newid += 1
            newtxt = "NEW{0:05d}".format(newid)
            tree_labels.append(tree_txt + newtxt) # Add the name and make a new label
            treeall_labels.append(tree_txt + newtxt)

            # Delete rows/columns/labels in "plain" distance matrix
            for i in sorted(min_ind, reverse=True):
                del tree_labels[i]
                newnojoin = []
                for pair in nojoin:
                    p1, p2 = pair
                    if p1 > min_ind[0]:
                        p1 -= 1
                    if p2 > min_ind[0]:
                        p2 -= 1
                    if p1 > min_ind[1]:
                        p1 -= 1
                    if p2 > min_ind[1]:
                        p2 -= 1
                    newnojoin.append((p1, p2))
                nojoin = newnojoin
            work_mx = np.delete(work_mx, min_ind, 0)  # Delete the 2 rows
            work_mx = np.delete(work_mx, min_ind, 1)  # Delete the 2 columns

            for i in sorted(min_ind_all, reverse=True):
                if treeall_labels[i] not in labels:
                    del treeall_labels[i]
                    workall_mx = np.delete(workall_mx, [i], 0)
                    workall_mx = np.delete(workall_mx, [i], 1)
            

        """
        print("CHECK DOUBLES")
        epis = set()
        for l in tree_labels:
            es = [m.start() for m in re.finditer('EPI_', l)]
            print("STARTING POINTS", es)
            for i in es:
                s = ""
                print("NOW", i, l[i:])
                for c in l[i:]:
                    if c != ":":
                        s += c
                    else:
                        break
                print("CHECK FOR STRING", s)
                if s in epis:
                    print("DUPLICATE", s)
                else:
                    epis.add(s)
        """

#    print("HERE??")
#    print(tree_labels)
#    print(len(tree_labels))

    """
    es = [m.start() for m in re.finditer('EPI_', "({0}:{1:.2f},{2}:{3:.2f})ROOT;\n".format(tree_labels[0], max(0.01,newd_1), tree_labels[1], max(0.01,newd_2)))]
    print("STARTING POINTS", es)
    for i in es:
        s = ""
        print("NOW", i, l[i:])
        for c in l[i:]:
            if c != ":":
                s += c
            else:
                break
        print("CHECK FOR STRING", s)
        if s in epis:
            print("DUPLICATE", s)
        else:
            epis.add(s)
    """

    with open("PROVA.nwk", "w") as pf:
        if len(tree_labels) == 2:
            pf.write("({0}:{1:.5f},{2}:{3:.5f})ROOT;\n".format(tree_labels[0], max(THR,newd_1), tree_labels[1], max(THR,newd_2)))
        else:
            pf.write("{0};\n".format(tree_labels[0]))

    if len(tree_labels) == 2:
        return "({0}:{1:.5f},{2}:{3:.5f})ROOT;\n".format(tree_labels[0], max(THR, newd_1), tree_labels[1], max(THR,newd_2))
    else:
        return "{0};\n".format(tree_labels[0])

def test():
    d = [[0, 5,12,13],
        [ 5, 0,13,14], 
        [12,13, 0,11],
        [13,14,11, 0]]
    nhx = tree.nj(DistanceMatrix(d, ['c','d','f','g']), result_constructor=str)
    t = Tree(nhx)
    print(nhx)
    t.set_outgroup(t.get_midpoint_outgroup())
    print(t.get_ascii(show_internal=True))
    print(t.write(format=1))

    print("\n\nMY ALGO\n")
#    internal_nj(np.array(d), ['c','d','f','g'])

    print("CON NODI INTERNI")
    d = [[0, 1, 3, 4, 4, 9,10],
        [ 1, 0, 2, 3, 5,10,11],
        [ 3, 2, 0, 5, 7,12,13],
        [ 4, 3, 5, 0, 8,13,14],
        [ 4, 5, 7, 8, 0, 5, 6],
        [ 9,10,12,13, 5, 0,11],
        [10,11,13,14, 6,11, 0]]
    l = ['a','b','c','d','e','f','g']
    print(internal_nj(np.array(d), l))

    d = [[0.0000, 0.0001, 0.0003, 0.0004, 0.0004, 0.0009, 0.0010],
        [ 0.0001, 0.0000, 0.0002, 0.0003, 0.0005, 0.0010, 0.0011],
        [ 0.0003, 0.0002, 0.0000, 0.0005, 0.0007, 0.0012, 0.0013],
        [ 0.0004, 0.0003, 0.0005, 0.0000, 0.0008, 0.0013, 0.0014],
        [ 0.0004, 0.0005, 0.0007, 0.0008, 0.0000, 0.0005, 0.0006],
        [ 0.0009, 0.0010, 0.0012, 0.0013, 0.0005, 0.0000, 0.0011],
        [ 0.0010, 0.0011, 0.0013, 0.0014, 0.0006, 0.0011, 0.0000]]
    l = ['a','b','c','d','e','f','g']
    print(internal_nj(np.array(d), l))



if __name__ == "__main__":
    test()
