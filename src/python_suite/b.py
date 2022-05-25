import sys
A = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

with open('blosumtmp.txt') as bf:
    il = -1
    for line in bf:
        il += 1
        if not il:
            continue
        f = line.split()
        for ii, fi in enumerate(f[1:]):
            if A[il-1] != A[ii]:
                print("{0}\t{1}\t{2:4d}".format(A[il-1], A[ii], -int(fi)))
            else:
                print("{0}\t{1}\t{2:4d}".format(A[il-1], A[ii], 0))
## WITH GAPS:
        print("{0}\t{1}\t{2}".format(A[il-1], '-', float('inf')))
        print("{0}\t{1}\t{2}".format(A[il-1], 'X', float('inf')))
for a in A:
    print("{0}\t{1}\t{2}".format('-', a, -10))
    print("{0}\t{1}\t{2}".format('X', a, -10))
print("-\t-\t0")
print("X\t-\t0")
print("-\tX\t0")
print("X\tX\t0")
