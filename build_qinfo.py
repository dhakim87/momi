import pandas as pd
import os
import subprocess


taxids = pd.read_csv("./ncbi_taxids.txt",
                     sep="\t")

unique_ids = set(taxids['ncbi_id'].tolist())
unique_ids = sorted(list(unique_ids))
print("Num Ids:", len(unique_ids))
for i in range(len(unique_ids)):
    print(i, "/", len(unique_ids))
    ncbi_id = unique_ids[i]
    print(ncbi_id)
    x = subprocess.call("esearch -db protein -query txid" + str(ncbi_id) + "[Organism] > ./proteomes/" + str(ncbi_id) + ".qinfo",
                        shell=True)
