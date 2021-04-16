import os;
import sys;
import pandas as pd;
from collections import defaultdict

MIN_SCORE = int(sys.argv[1])
print("BLOSUM62 THRESHOLD >=", MIN_SCORE)

taxids = pd.read_csv("./ncbi_taxids.txt", sep="\t")

ncbi_to_genome = defaultdict(list)

genome_ids = [str(x) for x in taxids['genome_id'].tolist()]
ncbi_ids = [str(x) for x in taxids['ncbi_id'].tolist()]

for i in range(len(ncbi_ids)):
    ncbi_to_genome[ncbi_ids[i]].append(genome_ids[i])

print("Num Ids:", len(ncbi_ids))

files_processed = 0
with open("./mbp_table" + str(MIN_SCORE) + ".csv", "w") as out_file:
    for filename in os.listdir("./netmhc"):
        if not filename.endswith(".mbp"):
            continue
        
        files_processed += 1
        if files_processed % 500 == 0:
            print(files_processed);
        
        ncbi_id = filename.split(".")[0]
        putative_mimics = []
        with open("./netmhc/" + filename) as f:
            # Skip header
            line = f.readline()
        
            # Start parsing
            line = f.readline()
            while line != "":
                ss = line.split(",")
                mbp_score = int(ss[0].strip())
                mimic = ss[1].strip()
                if mbp_score >= MIN_SCORE:
                    putative_mimics.append(mimic)
                else:
                    break;
                line = f.readline()
        
        putative_mimics = sorted(putative_mimics)
        for genome_id in ncbi_to_genome[ncbi_id]:
            out_file.write(ncbi_id + "," + genome_id + "," + ";".join(putative_mimics) + "\n")
