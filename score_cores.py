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


def print_all(rows):
    for row in rows:
        print(", ".join([str(r) for r in row]))


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

    max_sum = max(sums)
    sum_index = sums.index(max_sum)

    return max_sum, core_epitopes[sum_index]


# Strong and weak cores predicted by NetMHCIIpan v4.0
MBP_CORES_STRONG = [
    "VHFFKNIVT"
]
MBP_CORES_WEAK = [
    "LIRLFSRDA",
    "LDVMASQKR",
    "IGRFFGGDR",
    "TAHYGSLPQ",
    "FKNIVTPRT",
    "IFKLGGRDS"
]
MOG_CORES_STRONG = [
    "IRALVGDEV",
    "VHLYRNGKD"
]
MOG_CORES_WEAK = [
    "VGWYRPPFS",
    "IVPVLGPLV"
]
PLP1_CORES_STRONG = [
    "IETYFSKNY",
    "VLPWNAFPG",
    "IAAFVGAAA"
]
PLP1_CORES_WEAK = [
    "INVIHAFQY",
    "YVIYGTASF",
    "FFLYGALLL",
    "LLLAEGFYT",
    "FGDYKTTIC",
    "VYIYFNTWT",
    "YIYFNTWTT",
    "WNAFPGKVC",
    "FHLFIAAFV",
    "IAATYNFAV"
]

if len(sys.argv) >= 3 and sys.argv[2] == "DEBUG":
    ALL_CORES = MBP_CORES_STRONG + MBP_CORES_WEAK + \
                MOG_CORES_STRONG + MOG_CORES_WEAK + \
                PLP1_CORES_STRONG + PLP1_CORES_WEAK
    for core in ALL_CORES:
        print(core)
        print(aligned_score([core], core))
    exit(-1)

if sys.argv[1] == "MBP":
    CORE_EPITOPES = MBP_CORES_STRONG + MBP_CORES_WEAK
    SCORE_NAME = "MBP_SCORE"
elif sys.argv[1] == "MOG":
    CORE_EPITOPES = MOG_CORES_STRONG + MOG_CORES_WEAK
    SCORE_NAME = "MOG_SCORE"
elif sys.argv[1] == "PLP1":
    CORE_EPITOPES = PLP1_CORES_STRONG + PLP1_CORES_WEAK
    SCORE_NAME = "PLP1_SCORE"
else:
    print("Unknown target protein:", sys.argv[1])
    raise Exception("Unknown target protein: " + sys.argv[1])

putative_mimics = set([])
for line in sys.stdin:
    putative_mimic = line.strip()
    putative_mimics.add(putative_mimic)

output_rows = []
for mimic in putative_mimics:
    if len(mimic) != len(CORE_EPITOPES[0]):
        continue
    score, core_epitope = aligned_score(CORE_EPITOPES, mimic)
    output_rows.append((score, mimic, core_epitope))
output_rows.sort(key=lambda row: row[0], reverse=True)

print(SCORE_NAME + ", MIMIC, EPITOPE")
print_all(output_rows)
