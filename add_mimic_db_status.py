import sqlite3
import subprocess
import os
import xml.etree.ElementTree as ET
from os import path
import pandas as pd

def parse_qinfo(ncbi_id):
    tree = ET.parse("proteomes/" + ncbi_id + ".qinfo")
    root = tree.getroot()
    count = root.findall("Count")[0]
    return int(count.text)

def count_proteins_retrieved(ncbi_id):
    if not path.exists("proteomes/" + ncbi_id + ".fasta"):
        return -1
    result = subprocess.run("grep '>' proteomes/" + ncbi_id + ".fasta | wc -l",
                            stdout=subprocess.PIPE, shell=True)
    count = int(result.stdout.decode().strip())
    return count

taxids = pd.read_csv("./ncbi_taxids.txt", sep="\t")
ncbi_ids = set([str(x) for x in taxids['ncbi_id'].tolist()])
print("Num Ids:", len(ncbi_ids))

retrieved_ncbi_ids = set()
for file in os.listdir("proteomes"):
    if file.endswith(".qinfo"):
        retrieved_ncbi_ids.add(file.split('.')[0])

conn = sqlite3.connect("mimics.db")
c = conn.cursor()
c.execute("DROP TABLE status")
c.execute("CREATE TABLE IF NOT EXISTS status (ncbi_id text PRIMARY KEY, file_status TEXT, processing_status TEXT)")

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

    status = "valid"
    if actual_proteins == -1:
        status = "no_file"
    elif actual_proteins < expected_proteins:
        status = "missing_proteins"
    elif actual_proteins == expected_proteins:
        status = "valid"
    else:
        status = "wtf" # Garbage tools downloaded extra proteins...
    
    c.execute("INSERT INTO status(ncbi_id, file_status) VALUES(?,?)",(ncbi_id, status))




conn.commit()
conn.close()

