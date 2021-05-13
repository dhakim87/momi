import os;
from os import path
import subprocess;
import xml.etree.ElementTree as ET
import sys;
import pandas as pd;
import sqlite3;
import csv

from collections import defaultdict

conn = sqlite3.connect("mimics.db")
c = conn.cursor()

def make_tables():
    c.execute("CREATE TABLE mimic (ncbi_id text, genome_id text, epitope text, mimicked text, mimicked_epitope text, blosum real)")
    c.execute("CREATE INDEX mimic_any ON mimic(blosum, epitope)")
    c.execute("CREATE INDEX mimic_ncbi ON mimic(ncbi_id, genome_id, blosum)")
    c.execute("CREATE TABLE genome (genome_id text, ncbi_id text, genus text, species text)")
    c.execute("CREATE INDEX genome_idx ON genome(genome_id, genus, species)")
    c.execute("CREATE INDEX genome_ncbi ON genome(ncbi_id, genus, species)")
    c.execute("CREATE TABLE status (ncbi_id text, file_status text, processing_status text)")

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


MIN_SCORE = int(sys.argv[1])
print("BLOSUM62 THRESHOLD >=", MIN_SCORE)

taxids = pd.read_csv("./ncbi_taxids.txt", sep="\t")

ncbi_to_genome = defaultdict(list)
genome_to_ncbi = defaultdict(str)
genome_ids = [str(x) for x in taxids['genome_id'].tolist()]
ncbi_ids = [str(x) for x in taxids['ncbi_id'].tolist()]

for i in range(len(ncbi_ids)):
    ncbi_to_genome[ncbi_ids[i]].append(genome_ids[i])
    genome_to_ncbi[genome_ids[i]] = ncbi_ids[i]

print("Num Ids:", len(ncbi_ids))

def set_mimics():
    files_processed = 0
    for filename in os.listdir("./netmhc"):
        mimicked = None
        if filename.endswith(".mbp"):
            mimicked = "MBP"
        elif filename.endswith(".mog"):
            mimicked = "MOG"
        elif filename.endswith(".plp1"):
            mimicked = "PLP1"
        else:
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
                score = int(ss[0].strip())
                mimic = ss[1].strip()
                epitope = ss[2].strip()

                if score >= MIN_SCORE:
                    putative_mimics.append((mimic, score, epitope))
                else:
                    break;
                line = f.readline()

        putative_mimics = sorted(putative_mimics)
        for genome_id in ncbi_to_genome[ncbi_id]:
            for tuple in putative_mimics:
                c.execute("INSERT INTO mimic VALUES(?,?,?,?,?,?)", (ncbi_id, genome_id, tuple[0], mimicked, tuple[2], tuple[1]))

def set_genome():
    with open("./wol_metadata.tsv") as wol_meta:
        mreader = csv.reader(wol_meta, delimiter='\t', quotechar='|')
        headers = next(mreader, None)
        id = headers.index("#genome")
        si = headers.index("species")
        gi = headers.index("genus")
        for r in mreader:
            genome_id = r[id]
            ncbi_id = genome_to_ncbi[genome_id]
            species = r[si]
            genus = r[gi]
            c.execute("INSERT INTO genome VALUES(?,?,?,?)", (genome_id, ncbi_id, genus, species))

def set_status():
    # Loop over ids to see what we have
    valid = []
    no_file = []
    missing_proteins = []
    just_wtf = []
    percentages = []

    for i in range(len(ncbi_ids)):
        ncbi_id = ncbi_ids[i]
        print(i, '/', len(ncbi_ids), ncbi_id)
        expected_proteins = parse_qinfo(ncbi_id)
        actual_proteins = count_proteins_retrieved(ncbi_id)

        file_status = "valid"
        if actual_proteins == -1:
            file_status = "no_file"
        elif actual_proteins < expected_proteins:
            file_status = "missing_proteins"
        elif actual_proteins == expected_proteins:
            file_status = "valid"
        else:
            file_status = "wtf" # Garbage tools downloaded extra proteins...

        processing_status = []
        if path.exists("./netmhc/" + ncbi_id + ".fin"):
            processing_status.append("mhc")
        if path.exists("./netmhc/" + ncbi_id + ".mbp-scored"):
            processing_status.append("mbp")
        if path.exists("./netmhc/" + ncbi_id + ".mog-scored"):
            processing_status.append("mog")
        if path.exists("./netmhc/" + ncbi_id + ".plp1-scored"):
            processing_status.append("plp1")

        c.execute("INSERT INTO status(ncbi_id, file_status, processing_status) VALUES(?,?,?)",(ncbi_id, file_status, "|".join(processing_status)))

make_tables()
set_mimics()
set_genome()
set_status()

conn.commit()
conn.close()

