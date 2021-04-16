import subprocess
import os
import pandas as pd
import xml.etree.ElementTree as ET
import subprocess
from os import path

def parse_qinfo(ncbi_id):
    tree = ET.parse("proteomes/" + ncbi_id + ".qinfo")
    root = tree.getroot()
    count = root.findall("Count")[0]
    return int(count.text)


def count_proteins_retrieved(ncbi_id):
    if not path.exists("proteomes/" + ncbi_id + ".fasta"):
        return -1
    result = subprocess.run("grep '>' proteomes/" + ncbi_id + ".fasta | wc -l", stdout=subprocess.PIPE, shell=True)
    count = int(result.stdout.decode().strip())
    return count

# Retrieve list of required ncbi ids from the web of life lookup table
taxids = pd.read_csv("./ncbi_taxids.txt",
                     sep="\t")

ncbi_ids = set([str(x) for x in taxids['ncbi_id'].tolist()])
print("Num Ids:", len(ncbi_ids))

# Retrieve the list of .qinfo files we've managed to download
retrieved_ncbi_ids = set()
for file in os.listdir("proteomes"):
    if file.endswith(".qinfo"):
        retrieved_ncbi_ids.add(file.split('.')[0])

# Fail out if we still need to get more qinfo files
got_all_info = retrieved_ncbi_ids == ncbi_ids
print("Retrieved All NCBI ID info files:", got_all_info)
if not got_all_info:
    exit(0)

# Loop over ids to see what we have
valid = []
no_file = []
missing_proteins = []
just_wtf = []
percentages = []

ncbi_ids = sorted(list(ncbi_ids))
for i in range(len(ncbi_ids)):
    ncbi_id = ncbi_ids[i]
    print(i, '/', len(ncbi_ids), ncbi_id)
    expected_proteins = parse_qinfo(ncbi_id)
    actual_proteins = count_proteins_retrieved(ncbi_id)

    if actual_proteins == -1:
        no_file.append(ncbi_id)
    elif actual_proteins < expected_proteins:
        missing_proteins.append(ncbi_id)
    elif actual_proteins == expected_proteins:
        valid.append(ncbi_id)
    else:
        just_wtf.append(ncbi_id)

    if expected_proteins > 0 and actual_proteins >= 0:
        percentages.append(actual_proteins / expected_proteins)
    else:
        percentages.append(0.0)


percent_df = pd.DataFrame(percentages, index=ncbi_ids, columns=["percentage"])
percent_df.to_csv("./percentages.tsv", sep="\t")
print("Download percentages written to: percentages.tsv")

print("---")
print("Missing Proteins: ")
print(missing_proteins)

print("Corrupt Files: ")
print(just_wtf)

print("---")
print("Stats: ")
print("Valid:", str(len(valid)))
print("No File:", str(len(no_file)))
print("Corrupt (Missing Proteins):", str(len(missing_proteins)))
print("Corrupt (Extra Proteins!?!  wtf are these?)", str(len(just_wtf)))


to_download = sorted(no_file + missing_proteins + just_wtf)

df = pd.DataFrame(to_download, columns=["ncbi_id"])
df.to_csv("./remaining_ncbi_ids.tsv", sep="\t")
print("Downloads required written to: remaining_ncbi_ids.tsv")

