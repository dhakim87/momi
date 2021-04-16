import os;
import sys;
import pandas as pd;
import sqlite3;
import csv

from collections import defaultdict

conn = sqlite3.connect("mimics.db")
c = conn.cursor()
#c.execute("ALTER TABLE genome ADD COLUMN ncbi_id text")

taxids = pd.read_csv("./ncbi_taxids.txt", sep="\t")

ncbi_to_genome = defaultdict(list)
genome_to_ncbi = defaultdict(str)
genome_ids = [str(x) for x in taxids['genome_id'].tolist()]
ncbi_ids = [str(x) for x in taxids['ncbi_id'].tolist()]

for i in range(len(ncbi_ids)):
    ncbi_to_genome[ncbi_ids[i]].append(genome_ids[i])
    if genome_to_ncbi[genome_ids[i]]:
        print(genome_to_ncbi[genome_ids[i]])
        raise Exception("Genome maps to two ncbi ids!")
    genome_to_ncbi[genome_ids[i]] = ncbi_ids[i]
    c.execute("UPDATE genome SET ncbi_id=? WHERE genome_id=?", (ncbi_ids[i], genome_ids[i]))

conn.commit()
conn.close()

