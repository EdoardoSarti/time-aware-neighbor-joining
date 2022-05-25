import sys

countries = []
with open("clist.txt") as cf:
    for line in cf:
        countries.append(line.strip())

countries_in_aln = {}
nope = {}
with open(sys.argv[1]) as af:
    for line in af:
        if line.startswith(">"):
            country = line.split('|')[-1].strip()
            if country in countries:
                if country not in countries_in_aln:
                    countries_in_aln[country] = 0
                countries_in_aln[country] += 1
            else:
                if country not in nope:
                    nope[country] = 0
                nope[country] += 1

s = 0
for c in sorted(list(countries_in_aln)):
    print(c, countries_in_aln[c])
    s += countries_in_aln[c]
print("TOTAL", s)
print("")

s = 0
for c in sorted(list(nope), key=lambda x : nope[x]):
    print(c, nope[c])
    s += nope[c]
print("TOTAL", s)
