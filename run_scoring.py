import sys
import pandas as pd
import os
import subprocess
import xml.etree.ElementTree as ET

print("Script:", sys.argv[0])
print("Job Index:", sys.argv[1])
print("Num Jobs:", sys.argv[2])
print("Core Set (MBP/MOG/PLP1):", sys.argv[3]) 

JOB_INDEX = sys.argv[1]
NUM_JOBS = sys.argv[2]

# Retrieve list of required ncbi ids from the web of life lookup table
taxids = pd.read_csv("./ncbi_taxids.txt",
                     sep="\t")

ncbi_ids = set([str(x) for x in taxids['ncbi_id'].tolist()])
print("Num Ids:", len(ncbi_ids))

ncbi_ids = sorted(list(ncbi_ids))
for i in range(int(JOB_INDEX), len(ncbi_ids), int(NUM_JOBS)):
    ncbi_id = ncbi_ids[i]
    print(i, '/', len(ncbi_ids), ncbi_id)
    print(ncbi_id)

    if not os.path.exists("./netmhc/" + str(ncbi_id) + ".done"):
        print("Skipping " + str(ncbi_id) + " (No .done)")
        continue

    if sys.argv[3] == "MBP":
        ending = ".mbp"
    elif sys.argv[3] == "MOG":
        ending = ".mog"
    elif sys.argv[3] == "PLP1":
        ending = ".plp1"
    else:
        raise Exception("Unknown Protein: " + sys.argv[3])

    if os.path.exists("./netmhc/" + str(ncbi_id) + ending + "-scored"):
        print("Skipping " + str(ncbi_id) + " (Existence of .scored says it's already scored)")
        continue

    #TODO: Decide between multiple sorted scoring files (.mbp, .mog, .plp) and some combined file or combined score
    cmd = subprocess.run('cat ./netmhc/' + str(ncbi_id) + '.core | python score_cores.py ' + sys.argv[3] + ' > ./netmhc/' + str(ncbi_id) + ending,
                         shell=True)
    cmd = subprocess.run("touch ./netmhc/" + str(ncbi_id) + ending + "-scored", shell=True)
