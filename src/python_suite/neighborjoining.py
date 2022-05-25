from skbio import DistanceMatrix
from skbio import tree
from ete3 import Tree
import numpy as np
import re

def consensus_sequence(seq_1, seq_2):
    cons = []
    for ic in range(len(seq_1)):
        if seq_1[ic] == "X" or seq_2[ic] == "X":
            consc = "X"
        elif seq_1[ic] == seq_2[ic]:
            consc = seq_1[ic]
        elif seq_1[ic] == "-":
            consc = seq_2[ic]
        elif seq_2[ic] == "-":
            consc = seq_1[ic]
        else:
            consc = "X"
        cons.append(consc)
    return "".join(cons)


def internal_nj(dist_mx, labels, sequences):
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

    # Initialize connected component sets
    conncomps = [{x} for x in labels]
    cc_dict = {x : i for i, x in enumerate(labels)}

    # Create working distance matrix
    row_labels = labels[:] # Labels of the plain distance / NJ matrix
    col_labels = labels[:] # Labels of the workall matrix (leaves are kept)
    work_mx = dist_mx.copy()  # The distance matrix whose columns never get reduced 

    newid = -1
    PRECISION = 5
    THR = 0.1**PRECISION
    while work_mx.shape[0] > 2 and len(conncomps) > 1:
        n, m = work_mx.shape

        print("Connected components", conncomps)
##        print("CC DICT", cc_dict)

        # REVIEW
        # Creating the parent dictionary
        parent = {} # Parent of each column label of work_mx
        for label in labels:
            for i, l in enumerate(row_labels):
                if l == label or ((label != l.split(")")[-1] and "TMP"+label != l.split(")")[-1] and "GNJ"+label != l.split(")")[-1]) and (("("+label+":" in l) or (","+label+":" in l) or (")"+label+":" in l))):
#                    print(col_labels[i], "IS PARENT OF", label)
                    parent[label] = i
                    break
            if label not in parent:
                for i, l in enumerate(col_labels):
                    if l.split(")")[-1].replace("TMP","").replace("GNJ","")  == label:
                        parent[label] = i
#                        print(col_labels[i], "IS PARENT OF", label)
                        break
                

        # Generalized Neighbor Joining matrix (NJ plus rectangular complement) 
        si = np.sum(work_mx[:,:n], axis=1) # Sum on columns, leave rows as variable
        gnj_mx = 99999.*np.ones(work_mx.shape)
        for i in range(work_mx.shape[0]):
            ni = col_labels[i].split(")")[-1].replace("TMP", "").replace("GNJ", "")
            # Plain NJ matrix
            for j in range(work_mx.shape[0]):
                if i == j or (i in is_gnj and j in is_gnj):
                    continue
                nj = col_labels[j].split(")")[-1].replace("TMP", "").replace("GNJ", "")
                if cc_dict[ni] == cc_dict[nj]:  # Same connected component!
                    continue
                gnj_mx[i,j] = (n-2)*work_mx[i,j] - si[i] - si[j] 
            # Complement
            if i in is_gnj:
                continue
            for j in range(work_mx.shape[0], work_mx.shape[1]):
                nj = col_labels[j].split(")")[-1].replace("TMP", "").replace("GNJ", "")
                if cc_dict[ni] == cc_dict[nj]:  # Same connected component!
                    continue
                gnj_mx[i,j] = 2*(n-2)*work_mx[i,j] - 2*si[parent[col_labels[j]]] - 2*(n-2)*work_mx[i,parent[col_labels[j]]]

        print("GNJ MX")
        print(gnj_mx)

        # Find minimum of the current GNJ matrix
        min_ind = np.unravel_index(np.argmin(gnj_mx, axis=None), gnj_mx.shape)
        print("MOVE", min_ind,  col_labels[min_ind[0]], col_labels[min_ind[1]])

        # Update connected components
        left_name = col_labels[min_ind[0]].split(")")[-1].replace("TMP","").replace("GNJ", "")
        right_name = col_labels[min_ind[1]].split(")")[-1].replace("TMP","").replace("GNJ", "")
        left_status = col_labels[min_ind[0]].split(")")[-1][:3]
        if left_status not in ["TMP", "GNJ"]:
            left_status = ""
        right_status = col_labels[min_ind[1]].split(")")[-1][:3]
        if right_status not in ["TMP", "GNJ"]:
            right_status = "" 
        inds, targets, stati = [], [left_name, right_name], [left_status, right_status]
        for target in targets:
            for icc, cc in enumerate(conncomps):
                if target in cc:
                    inds.append(icc)
        newset = conncomps[inds[0]] | conncomps[inds[1]]
        del conncomps[max(inds)]
        del conncomps[min(inds)]
        conncomps.append(newset)
        cc_dict = {} # The real data str for CCs is this dictionary, the set is only for clarity
        for icc, cc in enumerate(conncomps):
            for l in cc:
                cc_dict[l] = icc
        
        # Is it a new move? 
        is_new_tmp, is_new_gnj = False, False
        search_for_accessible_leaves = False
        if min_ind[1] < gnj_mx.shape[0]:
            # Plain NJ criterion
            is_new_gnj = False
            newd_1 = 0.5*work_mx[min_ind] + 1/(2*(n-2))*(si[min_ind[0]] - si[min_ind[1]])
            newd_2 = work_mx[min_ind] - newd_1
            if ((newd_1 <= THR and newd_2 > 2*THR and row_labels[min_ind[0]] in labels) or (newd_2 <= THR and newd_1 > 2*THR and row_labels[min_ind[1]] in labels)):
                is_new_tmp = True
            else:
                search_for_accessible_leaves = True
        else:
### MODIFICA QUI QUANDO UNO DEI DUE E' GIA' TMP DEVE DIVENTARE (a, b)TMPb + c ==> (a, b)GNJb + (b, c)TMPc
            # Creation of a degeneracy
            is_new_gnj = True
            newd_1 = work_mx[min_ind]
            newd_2 = 0	# At distance 0 from the old leaf (now root of a new subtree)
        newd = np.array([newd_1, newd_2])


        # ARRIVATO QUI
        if is_new_gnj:
            # a + (leaf) b -> (a:N, b:0)GNJb
            parent_name = parent[rightname] 
            connector = rightname
            child = row_labels[min_ind[0]]
            new_label = "({0}:{1:.5f},{2}:{3:.5f})GNJ{2}".format(row_labels[min_ind[0]], newd_1, rightname, 0)
            is_leaf = True
        elif is_new_tmp:
            # Definitions: what is degenerate and what isn't depends on the distances
            if newd_1 <= THR and newd_2 > 2*THR:
                degenerate_id = 0
                nondegenerate_id = 1
            else:
                degenerate_id = 1
                nondegenerate_id = 0
            degenerate_name = targets[degenerate_id]
            nondegenerate_name = targets[nondegenerate_id]
            degenerate_status = stati[degenerate_id]
            nondegenerate_status = stati[nondegenerate_id]

            # Creation of the new label
            is_chain, is_resolution = False, False
            if degenerate_status or nondegenerate_status:
                # There is at least one TMP or GNJ in the two names
                t = Tree(row_labels[min_ind[nondegenerate_id]]+";", format=1)
                child1_name, child2_name = [x.name for x in t.children]
                if child1_name == nondegenerate_name:
                    ex_degenerate_name = child1_name
                    ex_nondegenerate_label = row_label[mind_ind[1]]
                    ex_dist = t.get_distance(t.children[0])
                elif child2_name == nondegenerate_name:
                    ex_degenerate_name = child2_name
                    ex_nondegenerate_label = row_label[mind_ind[0]]
                    ex_dist = t.get_distance(t.children[1])
                else:
                    print("DEGENERACY ERROR")
                    exit(1)

                if nondegenerate_status == "TMP":
                    # Chain: previous TMP becomes a GNJ because it now has an ancestor but it still lacks one child
                    # (a1, a2)TMPa2 + °b -> (a1, a2)GNJa2 + (a2, b)TMPb
                    new_label_gnj = "({0}:{1:.5f}, {2}:{3:.5f})GNJ{2}".format(ex_nondegenerate_label, ex_dist, ex_degenerate_name, 0)
                    new_label_tmp = "({0}:{1:.5f}, {2}:{3:.5f})TMP{2}".format(ex_degenerate_name, newd[nondegenerate_id], degenerate_name, 0)
                elif degenerate_status == "TMP":
                    # TMP resolution: one TMP gets completed by the new addition
                    # °(a1, a2)TMPa2 + b -> (a1, b)a2
                    new_label = "({0}:{1:.5f},{2}:{3:.5f}){4}".format(ex_nondegenerate_label, ex_dist, nondegenerate_name, newd[nondegenerate_id], degenerate_name)                
                elif nondegenerate_status == "GNJ" or degenerate_status == "GNJ":
                    # GNJ resolutions have priority over TMPs and are thus treated in the Resolution section
                    # (a1, a2)GNJa2 + °b -> resolution (a1, b)a2
                    # °(a1, a2)GNJa2 + b -> resolution (a1, b)a2
                    is_resolution = True
            else:
                # Names are "simple": a new TMP is created
                new_label = "({0}:{1:.5f}, {2}:{3:.5f})TMP{2}".format(row_labels[nondegenerate_id], newd[nondegenerate_id], degenerate_name, 0)
        if is_resolution:
            # GNJ resolution: GNJ entry gets completed and then added back to its ancestor
            # (a1, a2)GNJa2 + b, (a2, c... -> (a1, b)a2, (a2, c... -> ((a1, b)a2, c...
            if left_status == "GNJ":
                gnj_id = 0
                no_gnj_id = 1
            else:
                gnj_id = 1
                no_gnj_id = 0
            gnj_name = targets[gnj_id]
            t = Tree(row_labels[min_ind[gnj_id]]+";", format=1)
            child1_name, child2_name = [x.name for x in t.children]
            if child1_name == gnj_name:
                connector_name = child1_name
                terminal_label = row_labels[min_ind[1]]
                ex_dist = t.get_distance(t.children[1])
            elif child2_name == gnj_name:
                connector_name = child2_name
                terminal_label = row_labels[min_ind[0]]
                ex_dist = t.get_distance(t.children[0])
            else:
                print("DEGENERACY ERROR GNJ")
                exit(1)
            replace_label = "({0}:{1:.5f},{2}:{3:.5f}){4}".format(terminal_label, ex_dist, row_labels[min_ind[no_gnj]], newd[no_gnj], connector_name)
        if not (is_new_gnj or is_new_tmp or is_resolution):
            # a + b -> (a, b)X
            # (a1, a2)TMPa2 + (b1, b2)TMPb2 -> ((a1)a2, (b1)b2)X
            is_nj = True
            if stati[0] and stati[1]:
                # Two TMP make a plain NJ, but with minimal distances and creating two new GNJs (one per side)
                t1 = Tree(row_labels[min_ind[0]]+";", format=1)
                t2 = Tree(col_labels[min_ind[1]]+";", format=1)
                new_child1_name = "{0}:{1:.5f}".format(t1.name, 0.00001) # Distances set to zero for this kind of operation
                new_child2_name = "{0}:{1:.5f}".format(t2.name, 0.00001)
                if t1.get_distance(t1.children[0]) > t1.get_distance(t1.children[1]):
                    new_grandchild1_label = "{0}".format(t1.children[0].write(format=1)[:-1])
                    ex_dist1 = t1.get_distance(t1.children[0])
                else:
                    new_grandchild1_label = "{0}".format(t1.children[1].write(format=1)[:-1])
                    ex_dist1 = t1.get_distance(t1.children[1])
                if t2.get_distance(t2.children[0]) > t2.get_distance(t2.children[1]):
                    new_grandchild2_label = "{0}".format(t2.children[0].write(format=1)[:-1])
                    ex_dist2 = t2.get_distance(t2.children[0])
                else:
                    new_grandchild2_label = "{0}".format(t2.children[1].write(format=1)[:-1])
                    ex_dist2 = t2.get_distance(t2.children[1])
                new_label = "({0}:{1:.5f}, {2}:{3:.5f})".format(new_child1_name, 0.00001, new_child2_name, 0.00001)  # Join the two nodes 
                new_gnj1_label = "({0}:{1:.5f}, {2}:{3:.5f})GNJ{2}".format(new_grandchild1_label, ex_dist1, new_child1_name, 0)
                new_gnj2_label = "({0}:{1:.5f}, {2}:{3:.5f})GNJ{2}".format(new_grandchild2_label, ex_dist2, new_child2_name, 0)
            elif (not stati[0]) and (not stati[1]):
                # This is the plain NJ
                new_label = "({0}:{1:.5f}, {2}:{3:.5f})".format(row_labels[min_ind[0], newd[0], row_labels[min_ind[1], newd[1]])
            else:
                print("NJ ERROR")
                exit(1)
            ancestor_name =
            


        # Check if new node is degenerate and write new part of text to add to tree
        degenerate = True
        if new_move:
            tree_txt = "({0}:{1:.5f},{2}:{3:.5f})".format(row_labels[min_ind[0]], newd_1, col_labels[min_ind[1]], newd_2)
            degenerate = False
        elif row_labels[min_ind[0]].split(")")[-1][:3] == "TMP" and col_labels[min_ind[1]].split(")")[-1][:3] == "TMP": # Only two TMP can merge
            t1 = Tree(row_labels[min_ind[0]]+";", format=1)
            t2 = Tree(col_labels[min_ind[1]]+";", format=1)
            new_child1 = "{0}:{1:.5f}".format(t1.name, 0.00001) # Distances set to zero for this kind of operation
            new_child2 = "{0}:{1:.5f}".format(t2.name, 0.00001)
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
        elif row_labels[min_ind[0]].split(")")[-1][:3] == "TMP" or row_labels[min_ind[0]].split(")")[-1][:3] == "GNJ": # TMP-normal and GNJ-normal
            t1 = Tree(row_labels[min_ind[0]]+";", format=1)
            for c in t1.children:
                if t1.get_distance(c) <= THR:
                    new_parent = c.name
                else:
                    new_child1 = "{0}".format(c.write(format=1)[:-1])
                    new_child2 = "{0}:{1:.5f}".format(col_labels[min_ind[1]], newd_1+newd_2)
            wdim, ldim = 0, 1
            tree_txt = "({0},{1}){2}".format(new_child1, new_child2, new_parent)  # Join the two nodes
        elif col_labels[min_ind[1]].split(")")[-1][:3] == "TMP" or col_labels[min_ind[1]].split(")")[-1][:3] == "GNJ": # normal-TMP and normal-GNJ
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
        

        print("TREETEXT", tree_txt)
        print(Tree(tree_txt+";", format=1).get_ascii(show_internal=True))
        print("NEW DISTANCES", newd)

        # Reconstruction of the consensus chain
        lname_1 = col_labels[min_ind[0]].split(")")[-1].replace("TMP","").replace("GNJ", "")
        lname_2 = col_labels[min_ind[1]].split(")")[-1].replace("TMP","").replace("GNJ", "")
        seq_1 = sequences[lname_1]
        seq_2 = sequences[lname_2]
        print("SEQUENCES")
        print(seq_1)
        print(seq_2)
        cons = consensus_sequence(seq_1, seq_2)

        # Search a possible leaf to associate to the internal node
        assigned = None
        if search_for_accessible_leaves and not new_move:
            accessible_leaves = [x for x in labels if (x in col_labels) and (x not in conncomps[-1])]
            #accessible_leaves = [x for x in labels if ("("+x+":" not in tree_txt) and (","+x+":" not in tree_txt) and (")"+x+":" not in tree_txt) and (x in col_labels)]  # internal nodes can be filled only by a subset of leaves
##            print("ACCESSIBLE LEAVES", accessible_leaves) 
            # Is there any ORIGINAL LEAF NOT IN THE TWO SUBTREES that respects these distances? Check
            # Calculate all distances between ORIGINAL LEAVES MINUS LEAVES IN SUBTREES and the new node
            # Use this array in the next loop, instead of work_mx
            t = []
            for label in accessible_leaves:
                icol = col_labels.index(label)
                lname = label.split(")")[-1].replace("TMP","").replace("GNJ", "")
                seq = sequences[lname]
                score = 0
                norm = 0
                for ic in range(len(seq_1)):
                    if (cons[ic] == seq[ic] and cons[ic] != "-" and cons[ic] != "X") or (cons[ic] == "X"):
#                    if cons[ic] == seq[ic] or cons[ic] == "-" or seq[ic] == "-" or cons[ic] == "X" or seq[ic] == "X":
                        score += 1
                    if cons[ic] != "-" and seq[ic] != "-" and seq[ic] != "X":
                        norm += 1
                if not norm:
                    score = 0
                else:
                    score /= norm
#                d = np.array([work_mx[min_ind[0], icol], work_mx[min_ind[1], icol]]) # search_for_accessible_leaves tag guarantees min_ind[1] < work_mx.shape[0]
#                t.append((np.linalg.norm(d - newd)/np.linalg.norm(newd), icol))  # Distances from the rest of the tree and difference between branch lengths
##                print(score)
                t.append((score, icol))
            t = sorted(t, key= lambda x: x[0])  # The difference between the branch length must be minimized. Optimal: 0, Threshold: 1 (difference must be less than original branch length)
     
            # Stub for more complex criterion of choice
            if t and t[0][0] > 0.8: # ACTIVATED
                assigned = t[0]


        # Manage matrices and label lists according to one of the five possibilities:
        if new_move:
            # 0) New degeneracy (right side)
            print("Child-parent degeneracy (GNJ move)")

            # Prepare text to be added
            tree_txt = "({0}:{1:.5f},{2}:{3:.5f})".format(row_labels[min_ind[0]], newd_1, col_labels[min_ind[1]], newd_2)

            # Bring the old leaf back to the NJ columns (as first)
            tmplab = tree_txt + "GNJ{0}".format(col_labels[min_ind[1]])
            no_gnj_labels.add("GNJ{0}".format(col_labels[min_ind[1]]))
            del col_labels[min_ind[1]]
            col_labels = [tmplab] + col_labels
            tmpcol = work_mx[:,min_ind[1]].copy()
            work_mx = np.delete(work_mx, [min_ind[1]], 1)
            work_mx = np.column_stack((tmpcol, work_mx))

            # Make the corresp row appear (and increase min_ind[0] by one!)
            row_labels = [tmplab] + row_labels
            tmprow = np.append(tmpcol, [0.]*(m-n-1))
            tmprow = np.append([0.], tmprow)
            work_mx = np.vstack((tmprow, work_mx))
            min_ind = (min_ind[0]+1, min_ind[1])

            # If left was simple leaf, it will be conserved
            if col_labels[min_ind[0]] in labels:
                work_mx = np.column_stack((work_mx, work_mx[:,min_ind[0]]))
                work_mx[0,-1] = 0.
                col_labels = col_labels + [col_labels[min_ind[0]]]
            work_mx = np.delete(work_mx, [min_ind[0]], 0)
            work_mx = np.delete(work_mx, [min_ind[0]], 1)
            del row_labels[min_ind[0]]
            del col_labels[min_ind[0]]
        elif row_labels[min_ind[0]] in labels and (newd_1 <= THR and newd_2 > 2*THR):
            # 1) New degeneracy (left side)
            print("Creation left-degeneracy")

            if row_labels[min_ind[0]].split(")")[-1][:3] == "TMP":
                # Normal pair
                
            elif row_labels[min_ind[1]].split(")")[-1][:3] == "TMP":
                # Chain
            elif row_labels[min_ind[0]].split(")")[-1][:3] == "GNJ": 
                # Normal pair, solve degeneration 
            elif row_labels[min_ind[1]].split(")")[-1][:3] == "GNJ":
                # Impossible: WARNING message, approximation
            else:
                # Creation of a TMP

            # Rename the row/column (internal)
            newtxt = "TMP{0}".format(col_labels[min_ind[0]])
            col_labels[min_ind[0]] = tree_txt + newtxt
            row_labels[min_ind[0]] = tree_txt + newtxt

            # If the non-degenerate node is a leaf, copy the column to onto the last position, then delete row and column in any case
            if col_labels[min_ind[1]] in labels:
                work_mx = np.column_stack((work_mx, work_mx[:,min_ind[1]]))
                col_labels = col_labels + [col_labels[min_ind[1]]]
            work_mx = np.delete(work_mx, [min_ind[1]], 1)
            del col_labels[min_ind[1]]
            if min_ind[1] < work_mx.shape[0]:
                work_mx = np.delete(work_mx, [min_ind[1]], 0)
                del row_labels[min_ind[1]]
        elif col_labels[min_ind[1]] in labels and (newd_1 > 2*THR and newd_2 <= THR):
            # 2) New degeneracy (right side)
            print("Creation right-degeneracy")

            # Rename the row/column (internal)
            newtxt = "TMP{0}".format(col_labels[min_ind[1]])
            col_labels[min_ind[1]] = tree_txt + newtxt
            if min_ind[1] < work_mx.shape[0]:
                row_labels[min_ind[1]] = tree_txt + newtxt
            work_mx = np.delete(work_mx, [min_ind[0]], 0)

            # If the non-degenerate node is a leaf, copy the column to onto the last position, then delete row and column in any case
            if col_labels[min_ind[0]] in labels:
                work_mx = np.column_stack((work_mx, work_mx[:,min_ind[0]]))
                col_labels = col_labels + [col_labels[min_ind[0]]]
            work_mx = np.delete(work_mx, [min_ind[0]], 1)
            del col_labels[min_ind[0]]
        elif row_labels[min_ind[0]].split(")")[-1][:3] == "TMP" or row_labels[min_ind[0]].split(")")[-1][:3] == "GNJ" or col_labels[min_ind[1]].split(")")[-1][:3] == "TMP" or col_labels[min_ind[1]].split(")")[-1][:3] == "GNJ":
            # 3) Previously declared degeneracy to be resolved
            print("Resolution of one degeneracy")

            # Define winner/loser
            if row_labels[min_ind[0]].split(")")[-1][:3] == "TMP" or row_labels[min_ind[0]].split(")")[-1][:3] == "GNJ":
                wdim, ldim = 0, 1
            else:
                wdim, ldim = 1, 0

            # Prepare text to be added
            t1 = Tree(row_labels[wdim]+";", format=1)
            for c in t1.children:
                if t1.get_distance(c) <= THR:
                    new_parent = c.name
                else:
                    new_child1 = "{0}".format(c.write(format=1)[:-1])
                    new_child2 = "{0}:{1:.5f}".format(col_labels[ldim], newd_1+newd_2)
            tree_txt = "({0},{1}){2}".format(new_child1, new_child2, new_parent)  # Join the two nodes

            # Replace in winner label, replace winner name with loser label
            llabel = col_labels[min_ind[ldim]]
            wname = col_labels[min_ind[wdim]].split(")")[-1][3:] # Record in col_labels has to end with ")TMPxxxx"
            col_labels[min_ind[wdim]] = (")".join(col_labels[min_ind[wdim]].split(")")[:-1])+")").replace(wname, llabel) + wname
            if min_ind[wdim] < work_mx.shape[0]:
                row_labels[min_ind[wdim]] = col_labels[min_ind[wdim]]

            # Delete loser, if it was a leaf then record it as last column
            if col_labels[min_ind[ldim]] in labels:
                work_mx = np.column_stack((work_mx, work_mx[:,min_ind[ldim]]))
                col_labels = col_labels + [col_labels[min_ind[ldim]]]
            work_mx = np.delete(work_mx, [min_ind[ldim]], 0)
            work_mx = np.delete(work_mx, [min_ind[ldim]], 1)
            del row_labels[min_ind[ldim]]
            del col_labels[min_ind[ldim]]
            if min_ind[ldim] < parent[wname]:
                parent[wname] -= 1

            # Replace parent column name
            col_names = [x.split(")")[-1].replace("TMP", "").replace("GNJ", "") for x in col_labels]
            newcol_idx = col_names.index(wname)
##            print("AAAAAA", parent[wname], newcol_idx, wname)
##            print(col_labels)
##            print(row_labels)
            if parent[wname] != newcol_idx:
                col_labels[parent[wname]] = col_labels[parent[wname]].replace(wname, col_labels[newcol_idx])
                row_labels[parent[wname]] = row_labels[parent[wname]].replace(wname, row_labels[newcol_idx])
                del col_labels[newcol_idx]
                del row_labels[newcol_idx]
                work_mx = np.delete(work_mx, [newcol_idx], 0)
                work_mx = np.delete(work_mx, [newcol_idx], 1)
        elif type(assigned) != type(None): 
            # 4) Leaf assigned to be the new internal node
            print("Internal node assignment")

            # No row/column to add: just one to rename and one to delete
            clave = col_labels[assigned[1]]
##            print("BEF", col_labels)
#            col_labels[assigned[1]] = tree_txt + col_labels[assigned[1]]
            if assigned[1] < work_mx.shape[0]: # If leaf is still to be positioned, we have to update rows too
                 row_labels[assigned[1]] = tree_txt + row_labels[assigned[1]]
                 col_labels[assigned[1]] = tree_txt + col_labels[assigned[1]]
            else:
                 for il, l in enumerate(row_labels):
                     if clave in l and l.split(")")[-1].replace("TMP", "").replace("GNJ", "") != clave:
                         break
                 row_labels[il] = row_labels[il].replace(clave, tree_txt + clave)
                 col_labels[il] = col_labels[il].replace(clave, tree_txt + clave)
##            print("AFT", col_labels)

            # Place columns last if they are leaves
            if row_labels[min_ind[0]] in labels:
                work_mx = np.column_stack((work_mx, work_mx[:,min_ind[0]]))
                col_labels = col_labels + [col_labels[min_ind[0]]]
            if col_labels[min_ind[1]] in labels:
                work_mx = np.column_stack((work_mx, work_mx[:,min_ind[1]]))
                col_labels = col_labels + [col_labels[min_ind[1]]]
            # Delete original rows and columns
            if min_ind[0] < min_ind[1]:
                del col_labels[min_ind[1]]
                del col_labels[min_ind[0]]
                del row_labels[min_ind[1]]
                del row_labels[min_ind[0]]
            else:
                del col_labels[min_ind[0]]
                del col_labels[min_ind[1]]
                del row_labels[min_ind[0]]
                del row_labels[min_ind[1]]
            work_mx = np.delete(work_mx, min_ind, 0)
            work_mx = np.delete(work_mx, min_ind, 1)
##            print("AFTAFT", col_labels)
            # Connected component update (new one is always last)
            # Merge connected component from the internal node
            for icc, cc in enumerate(conncomps):
                if clave in cc:
                    break
            conncomps[-1] |= conncomps[icc]
            del conncomps[icc]
            # Update all cc_dict after deleting one conncomp (shift in indexing)
            cc_dict = {}
            for icc, cc in enumerate(conncomps):
                for l in cc:
                    cc_dict[l] = icc


            # Delete assigned node from col_labels (not free anymore!)
            if clave in col_labels:
                icl = col_labels.index(clave)
                del col_labels[icl]
                work_mx = np.delete(work_mx, icl, 1)

        else:
            # 5) Classical neighbor joining
            print("Creation of new internal node (classic NJ)")
            # Plain NJ option
            # - Adds one row/column in 'plain' distances, associated to a new node, with the classical NJ formula
            # - Deletes the two rows/columns in 'plain' distances, associated to the two joined subtrees
            # - Adds one corresponding row/column in 'second' distance matrix, with the classical NJ formula (yet calculated on the 'second' matrix instead)

            # Distances new row/column
            u_col = np.zeros(work_mx.shape[0]+1)
            for i in range(1, work_mx.shape[0]+1):
                u_col[i] = 0.5*(work_mx[min_ind[0], i-1] + work_mx[i-1, min_ind[1]] - work_mx[min_ind[0], min_ind[1]]) # Classical NJ formula

            # Distances new row/column in second matrix
            u_row = np.zeros(work_mx.shape[1])
            u_row[:work_mx.shape[0]] = u_col[1:]
            for i in range(work_mx.shape[0], work_mx.shape[1]):
                name = col_labels[i].split(")")[-1].replace("TMP", "").replace("GNJ", "")
                # dist(new, oldleaf) = dist(new,parent(oldleaf)) + dist(parent(oldleaf), oldleaf) = NJdist(new,parent(oldleaf)) + dist(parent(oldleaf), oldleaf)
                u_row[i] = 0.5*(work_mx[min_ind[0], parent[name]] + work_mx[parent[name], min_ind[1]] - work_mx[min_ind[0], min_ind[1]]) + work_mx[parent[name], i]

            # Add row/column to matrices
            work_mx = np.vstack((u_row, work_mx))
            work_mx = np.column_stack((u_col, work_mx))

            # Update labels of matrices
            newid += 1
            newtxt = "NEW{0:05d}".format(newid)
            row_labels = [tree_txt + newtxt] + row_labels # Add the name and make a new label
            col_labels = [tree_txt + newtxt] + col_labels
            min_ind = (min_ind[0]+1, min_ind[1]+1)

            # New sequence
            sequences[newtxt] = cons

            # Delete rows, and place columns last if they are leaves, and delete them from where they are in any case 
            if row_labels[min_ind[0]] in labels:
                work_mx = np.column_stack((work_mx, work_mx[:,min_ind[0]]))
                col_labels = col_labels + [col_labels[min_ind[0]]]
            if col_labels[min_ind[1]] in labels:
                work_mx = np.column_stack((work_mx, work_mx[:,min_ind[1]]))
                col_labels = col_labels + [col_labels[min_ind[1]]]
            if min_ind[0] < min_ind[1]:
                del col_labels[min_ind[1]]
                del col_labels[min_ind[0]]
                del row_labels[min_ind[1]]
                del row_labels[min_ind[0]]
            else:
                del col_labels[min_ind[0]]
                del col_labels[min_ind[1]]
                del row_labels[min_ind[0]]
                del row_labels[min_ind[1]]
            work_mx = np.delete(work_mx, min_ind, 0)
            work_mx = np.delete(work_mx, min_ind, 1)

            # Update connected components
            conncomps[-1].add(newtxt)
            cc_dict[newtxt] = len(conncomps)-1

#    print("FINAL CCS", conncomps)
#    print("FINAL LABELS", row_labels)

    # Solve all remaining TMP and GNJ
    go = False
    for l in row_labels:
        if "TMP" in l.split(')')[-1] or "GNJ" in l.split(')')[-1]:
            go = True
            break
    while go:
#        print("SOLVE")
        for il, l in enumerate(row_labels):
            if "TMP" in l.split(')')[-1] or "GNJ" in l.split(')')[-1]:
                name = l.split(')')[-1].replace("GNJ", "").replace("TMP", "")
                txt = ",".join(l.split(',')[:-1])+")"+name
                break
#        print("NAME", name)
#        print("REPL TXT", txt)
        for il2, l in enumerate(row_labels):
            if name in l and name != l.split(')')[-1].replace("GNJ", "").replace("TMP", ""):
                break
#        print("LAB TO BE REPL", row_labels[il2])
        row_labels[il2] = row_labels[il2].replace(name, txt)
#        print("NEW LAB", row_labels[il2])
        del row_labels[il]
#        print(row_labels)
        go = False
        for l in row_labels:
            if "TMP" in l.split(')')[-1] or "GNJ" in l.split(')')[-1]:
                go = True
                break

    if len(row_labels) > 1:
        row_labels[0] = "({0}:{1:.5f}, {2}:{3:.5f})".format(row_labels[0], newd_1, row_labels[1], newd_2)
        del row_labels[1]

##    print("FINAL ROWS", row_labels)
#    print("FINAL COLS", col_labels)
#    print("FINAL TREE", tree_txt)
#    name = row_labels[0].split(")")[-1]
#    row_labels[0] = ")".join(row_labels[0].split(")")[:-1]) + "," + name + ":0.00001)"
    t = Tree("{0};\n".format(row_labels[0]), format=1)
    t.unroot()
##    print(t.get_ascii(show_internal=True))
    n = t.search_nodes(name="EPI_ISL_402124")[0]
#    print(n.get_ascii(show_internal=True))
#    print(t.get_ascii(show_internal=True))
    t.set_outgroup(n)
    print(t.get_ascii(show_internal=True))
    print("{0}".format(t.write(format=1)))
    with open("leg.txt", "w") as lf:
        lf.write("TREE_COLORS\nSEPARATOR TAB\nDATA\n") 
        for l in labels:
            if l[-5:] == "_THIS":
                lf.write("{0}\tbranch\trgb(204, 0, 0)\tnormal\t4\n".format(l))
    exit(1)

    # Last instance
    if len(row_labels) > 1 and len(conncomps) > 1:
        if row_labels[0].split(")")[-1][:3] == "TMP" and col_labels[1].split(")")[-1][:3] == "TMP":
            t1 = Tree(row_labels[0]+";", format=1)
            t2 = Tree(col_labels[1]+";", format=1)
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
            tree_txt = "(({0}){1}, ({2}){3})ROOT".format(new_grandchild1, new_child1, new_grandchild2, new_child2)  # Join the two nodes 
            degenerate = False
        elif row_labels[0].split(")")[-1][:3] == "TMP":
            t1 = Tree(row_labels[0]+";", format=1)
            for c in t1.children:
                if t1.get_distance(c) <= THR:
                    new_parent = c.name
                else:
                    new_child1 = "{0}".format(c.write(format=1)[:-1])
                    new_child2 = "{0}:{1:.5f}".format(col_labels[1], newd_1+newd_2)
            wdim, ldim = 0, 1
            tree_txt = "({0},{1}){2}".format(new_child1, new_child2, new_parent)  # Join the two nodes
        elif col_labels[1].split(")")[-1][:3] == "TMP":
            t1 = Tree(col_labels[1]+";", format=1)
            for c in t1.children:
                if t1.get_distance(c) <= THR:
                    new_parent = c.name
                else:
                    new_child1 = "{0}:{1:.5f}".format(row_labels[0], newd_1+newd_2)
                    new_child2 = "{0}".format(c.write(format=1)[:-1])
            wdim, ldim = 1, 0
            tree_txt = "({0},{1}){2}".format(new_child1, new_child2, new_parent)  # Join the two nodes
        else:
            degenerate = False
            new_parent = ""
            new_child1 = "{0}:{1:.5f}".format(row_labels[0], max(THR, newd_1))
            new_child2 = "{0}:{1:.5f}".format(col_labels[1], max(THR, newd_2))
            tree_txt = "({0},{1})ROOT".format(new_child1, new_child2, new_parent)  # Join the two nodes
    else:
        tree_txt = row_labels[0]
            
    with open("PROVA.nwk", "w") as pf:
        pf.write("{0};\n".format(tree_txt))

    print(Tree("{0};\n".format(tree_txt), format=1).get_ascii(show_internal=True))
    return "{0};\n".format(tree_txt)


def test():
    d = [[0, 1, 2, 3, 4, 5, 6],
        [ 1, 0, 1, 2, 3, 4, 5],
        [ 2, 1, 0, 1, 2, 3, 4],
        [ 3, 2, 1, 0, 1, 2, 3],
        [ 4, 3, 2, 1, 0, 1, 2],
        [ 5, 4, 3, 2, 1, 0, 1],
        [ 6, 5, 4, 3, 2, 1, 0]]    
    seqs = {
        'a' : 'AAAAAA',
        'b' : 'AAAAAC',
        'c' : 'AAAACC',
        'd' : 'AAACCC',
        'e' : 'AACCCC',
        'f' : 'ACCCCC',
        'g' : 'CCCCCC'
    }
    internal_nj(np.array(d), ['a','b','c','d','e','f','g'], seqs)

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
