import pandas as pd
import os
import subprocess
import xml.etree.ElementTree as ET


def parse_qinfo(ncbi_id):
    tree = ET.parse("proteomes/" + str(ncbi_id) + ".qinfo")
    root = tree.getroot()
    count = root.findall("Count")[0]
    return int(count.text)


def count_proteins_retrieved(ncbi_id):
    # Count proteins in a fasta file, 0 if file empty, -1 if file doesn't exist
    if not os.path.exists("proteomes/" + str(ncbi_id) + ".fasta"):
        return -1
    result = subprocess.run("grep '>' proteomes/" + str(ncbi_id) + ".fasta | wc -l", stdout=subprocess.PIPE, shell=True)
    count = int(result.stdout.decode().strip())
    return count


taxids = pd.read_csv("./ncbi_taxids.txt",
                     sep="\t")

# taxids = pd.read_csv("./remaining_ncbi_ids.tsv",
#                      sep="\t")

print("Your Path Is:")
x = subprocess.call("echo $PATH", shell=True)
print("Your API Key Is:")
x = subprocess.call("echo $NCBI_API_KEY", shell=True)

for ncbi_id in taxids["ncbi_id"]:
    print(ncbi_id)

    count = parse_qinfo(ncbi_id)
    # if count > 100000:
    #     print("Skipping " + str(ncbi_id) + " Too many proteins: (" + str(count) + ")")
    #     continue

    received = count_proteins_retrieved(ncbi_id)
    # if received > count * .5:
    #     print("Skipping " + str(ncbi_id) + " We already got more than half.  Uugghh")
    #     continue

    if received == count:
        print("Skip " + str(ncbi_id) + " File already downloaded")
        continue

    if received == 0:
        print("Retrying " + str(ncbi_id) + " (Which failed instantly in the past)")

    x = subprocess.call("esearch -db protein -query txid" + str(ncbi_id) + "[Organism] | efetch -format fasta > ./proteomes/" + str(ncbi_id) + ".fasta",
                        shell=True)
