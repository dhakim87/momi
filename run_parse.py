import sys
import pandas as pd
import os
import subprocess
import xml.etree.ElementTree as ET

print("Script:", sys.argv[0])
print("Job Index:", sys.argv[1])
print("Num Jobs:", sys.argv[2])

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

    if not os.path.exists("./netmhc/" + str(ncbi_id) + ".fin"):
        print("Skipping " + str(ncbi_id) + " (No .fin)")
        continue

    if os.path.exists("./netmhc/" + str(ncbi_id) + ".done"):
        print("Skipping " + str(ncbi_id) + " (Existence of .done says it's already parsed)")
        continue

    cmd = subprocess.run('grep "<=SB" ./netmhc/' + str(ncbi_id) + '.mhc | python parse_strong_binders.py 2> ./netmhc/' + str(ncbi_id) + '.parseErr > ./netmhc/' + str(ncbi_id) + '.core',
                           shell=True)
    cmd = subprocess.run("touch ./netmhc/" + str(ncbi_id) + ".done", shell=True)

