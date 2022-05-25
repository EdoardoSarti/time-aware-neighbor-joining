import sys

with open(sys.argv[1]) as f:
    iline = 0
    for line in f:
        iline += 1
        if iline < 4:
            print(line.strip())
        else:
            fields = line.split()
            print("\t".join(fields[:3]) + " " + " ".join(fields[3:5])[:-1] + ")\tnormal\t4")
