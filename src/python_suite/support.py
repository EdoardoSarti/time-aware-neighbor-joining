import os
import sys
from skbio import DistanceMatrix
from skbio import tree
import argparse
import numpy as np
import copy
from ete3 import Tree
import string
import random
import re

INF = float('Inf')

def errprint(text):
    print(text, file=sys.stderr)


def seqid(alnseq1, alnseq2, mode='symm-gap', keeplower=False):
    """
    symm: # matches / (# matches + # mismatches)
    symm-gap: # matches / min(length(seq))
    asymm: # matches / # length 1st seq
    """

    L = len(alnseq1)
    if len(alnseq2) != L:
        errprint("[seqid] FATAL: the two input sequences are not equally long")
        errprint("alnseq1")
        errprint(alnseq1)
        errprint("alnseq2")
        errprint(alnseq2)
        errprint("")
        exit(1)
    if alnseq1 != alnseq1.upper() or alnseq2 != alnseq1.upper():
        errprint("[seqid] WARNING: there are lowercase characters in the sequences.")
        errprint("keeplower option is {0}: such positions will {1}be disregarded".format("keeplower", "not "*(keeplower)))
        if keeplower:
            alnseq1 = alnseq1.replace('.', '-').upper()
            alnseq2 = alnseq2.replace('.', '-').upper()
    # TO DO: put a control on other symbols than [A-Za-z-.] 
    
    match, mismatch = 0, 0
    for i in range(L):
        if alnseq1[i] == alnseq2[i]:
            match += 1
        if alnseq1[i] != alnseq2[i] and alnseq1[i] not in ['.', '-'] and alnseq2[i] not in ['.', '-']:
            mismatch += 1
    if mode == 'symm':
        if match + mismatch == 0:
            return 0
        else:
            return match / (match + mismatch)
    elif mode == 'symm-gap':
        return match / min(len(alnseq1.replace('.','').replace('-','')), len(alnseq2.replace('.','').replace('-','')))
    elif mode == 'asymm':
        return match / len(alnseq1.replace('.','').replace('-',''))
    else:
        errprint("[seqid] FATAL: invalid mode: {0}".format(mode))
        exit(1)


def import_aln(aln_fn, aln_format=None):
    with open(aln_fn) as af:
        for line in af:
            if line.startswith(">") or aln_format.lower()=='fasta' or aln_format.lower()=='afa':
                return parser_fasta(aln_fn)
            elif line.startswith("#") or aln_format.lower()=='stockholm':
                return parser_stockholm(aln_fn)
            else:
                errprint("[import aln] FATAL: alignment format was not recognized.")
                errprint("recognized formats: fasta (afa), stockholm")


def parser_fasta(aln_fn, out_format="numpy"):
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

    if out_format == 'dict':
        return aln_seqs, names

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


def parser_stockholm(aln_fn):
    errprint("[parser_stockholm] FATAL: Not implemented yet")
    exit(1)


def adapt_aln_to_ref(aln, names, ref=0):
    M, L = aln.shape
    if M != len(names):
        errprint("[adapt_aln_to_ref] FATAL: alignment depth and number of labels are not the same")
        errprint("depth: {0}, len(names): {1}".format(M, len(names)))
        exit(1)
        
    if ref > len(names)-1:
        errprint("[adapt_aln_to_ref] FATAL: ref must be set in the interval [0, len(names)-1]")
        errprint("ref:", ref)
        errprint("len(names)", len(names))
        exit(1)

    del_columns = [i for i in range(L) if aln[ref][i] == '-' or aln[ref][i] == '.']
    del_aln = np.delete(aln, del_columns, 1)
    swap_names = names
    if ref != 0:
        del_aln[[0, ref]] = del_aln[[ref, 0]]
#        print(ref, len(names), names)
        swap_names = [names[ref]] + [names[i] for i in range(M) if i != ref]
    return del_aln, swap_names


def first_vs_all_seqid(aln):
    M, L = aln.shape
    return np.array([np.sum(aln[0,:] == aln[i,:])/L for i in range(M)])
