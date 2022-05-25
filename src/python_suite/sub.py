import sys

alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'j', 'k', 'i', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'v', 'w', 'x', 'y', 'z']

i = -1
with open(sys.argv[1]) as af:
    for line in af:
        if line.startswith('>'):
            i += 1
            fields = line.strip().split("|")
            fields[3] = alphabet[i//len(alphabet)//len(alphabet)] + alphabet[i//len(alphabet)%len(alphabet)] + alphabet[i%len(alphabet)]
#            fields[3] = alphabet[i%len(alphabet)]
            print("|".join(fields))
        else:
            print(line.strip())
