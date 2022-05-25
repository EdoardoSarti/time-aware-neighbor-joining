import sys

aln_fn = sys.argv[1]
list_fn = sys.argv[2]

if '--negative' in sys.argv:
    T, F = False, True
else:
    T, F = True, False

if '--onlynames' in sys.argv:
    onlynames = True
else:
    onlynames = False

names = []
with open(list_fn) as lf:
    for line in lf:
        names.append(line.split('|')[3])


with open(aln_fn) as af:
    go = False
    lastname = ""
    for line in af:
        if line.startswith(">"):
            if line.split('|')[3] in names:
                go = T
            else:
                go = F
            lastname = line.strip()
            continue
        if go:
            print(lastname)
            if not onlynames:
                print(line.strip())
            go = False
