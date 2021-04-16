import sys
import pandas as pd
import os
import subprocess
import xml.etree.ElementTree as ET

print("Script:", sys.argv[0])
print("Job Index:", sys.argv[1])
print("Num Jobs:", sys.argv[2])
print("NetMHCIIPan Location:", sys.argv[3])

JOB_INDEX = sys.argv[1]
NUM_JOBS = sys.argv[2]
NETMHCIIPAN = sys.argv[3]

if not os.path.exists(NETMHCIIPAN):
    print("Could not find NetMHCIIpan executable")
    exit(-1)

def parse_qinfo(ncbi_id):
    tree = ET.parse("proteomes/" + str(ncbi_id) + ".qinfo")
    root = tree.getroot()
    count = root.findall("Count")[0]
    return int(count.text)


def count_proteins_retrieved(ncbi_id):
    if not os.path.exists("proteomes/" + str(ncbi_id) + ".fasta"):
        return -1
    result = subprocess.run("grep '>' proteomes/" + str(ncbi_id) + ".fasta | wc -l", stdout=subprocess.PIPE, shell=True)
    count = int(result.stdout.decode().strip())
    return count


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

    count = parse_qinfo(ncbi_id)
#    if count > 50000:
#        print("Skipping " + str(ncbi_id) + " Too many proteins: (" + str(count) + ")")
#        continue

    if count == 0:
        print("Skipping " + str(ncbi_id) + " Empty File")
        continue

    if os.path.exists("./netmhc/" + str(ncbi_id) + ".fin"):
        print("Skipping " + str(ncbi_id) + " This has been previously attempted")
        continue

    received = count_proteins_retrieved(ncbi_id)

    if received < count * .95:
        print("Skipping " + str(ncbi_id) + " Less than 95% of proteins downloaded")
        continue

    cmd = subprocess.run(NETMHCIIPAN + " -filter -a DRB1_1501 -f ./proteomes/" + str(ncbi_id) + ".fasta > ./netmhc/" + str(ncbi_id) + ".mhc 2> ./netmhc/" + str(ncbi_id) + ".err",
                           shell=True)
    cmd = subprocess.run("touch ./netmhc/" + str(ncbi_id) + ".fin", shell=True)

