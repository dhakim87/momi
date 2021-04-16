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


def aligned_score(core_epitopes, putative_mimic):
    sums = []
    for core_epitope in core_epitopes:
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
        sums.append(sum)
    return max(sums)

MBP_CORES = ["VHFFKNIVT"]
MOG_CORES = ["IRALVGDEV", "VHLYRNGKD"]
PLP1_CORES = ["IETYFSKNY", "VLPWNAFPG", "IAAFVGAAA"]

if len(sys.argv) >= 3 and sys.argv[2] == "DEBUG":
    for core in MBP_CORES + MOG_CORES + PLP1_CORES:
        print(core)
        print(aligned_score([core],core))

if sys.argv[1] == "MBP":
    CORE_EPITOPES = MBP_CORES
    SCORE_NAME = "MBP_SCORE"
if sys.argv[1] == "MOG":
    CORE_EPITOPES = MOG_CORES
    SCORE_NAME = "MOG_SCORE"
if sys.argv[1] == "PLP1":
    CORE_EPITOPES = PLP1_CORES
    SCORE_NAME = "PLP1_SCORE"


cores = set([])
for line in sys.stdin:
    core = line.strip()
    cores.add(core)

mimics = [(aligned_score(CORE_EPITOPES, mimic), mimic) for mimic in cores if len(mimic) == len(CORE_EPITOPES[0])]
mimics = sorted(mimics, key=lambda mimic: mimic[0], reverse=True)

print(SCORE_NAME + ", CORE")
print_all(mimics)

