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
elif sys.argv[1] == "MOG":
    CORE_EPITOPES = MOG_CORES
    SCORE_NAME = "MOG_SCORE"
elif sys.argv[1] == "PLP1":
    CORE_EPITOPES = PLP1_CORES
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
