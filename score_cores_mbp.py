import sys

def readBlosum62():
    dict = {}
    with open("blosum62.txt") as f:
        data = f.read().splitlines()
        aas = data[0].split()

        i = 0
        for row in data[1:]:
            vals = row.split()[1:]
            for i in range(len(vals)):
                score = int(vals[i])
                colAA = aas[i]
                rowAA = row.split()[0]
                dict[(rowAA, colAA)] = score
        return dict

BLOSUM62 = readBlosum62()

def print_all(mimics):
    for mimic in mimics:
        row = str(mimic[0]) + "," + mimic[1]
        print(row)


def aligned_score(core_epitope, putative_mimic):
    sum = 0
    if len(core_epitope) != len(putative_mimic):
        raise Exception("Epitope length mismatch")
    for i in range(len(core_epitope)):
        key = (core_epitope[i], putative_mimic[i])
        if key in BLOSUM62:
            sum += BLOSUM62[key]
        else:
            # print("Unknown Amino Acid Match: ", key)
            sum -= 10
    return sum

MBP_CORE = "VHFFKNIVT"


cores = set([])
for line in sys.stdin:
    core = line.strip()
    cores.add(core)

mimics = [(aligned_score(MBP_CORE, mimic), mimic) for mimic in cores if len(mimic) == len(MBP_CORE)]
mimics = sorted(mimics, key=lambda mimic: mimic[0], reverse=True)

print("MBP_SCORE, CORE")
print_all(mimics)

