import json

from file_util import fasta_scan, fasta_glob

print("CONFIG INFO")
with open("run_config.json") as config_file:
    config = json.load(config_file)
    print(json.dumps(config, indent=4, sort_keys=True))


if "epitope_override" in config and len(config["epitope_override"]) > 0:
    print("TARGET EPITOPES:")
    for ep in config["epitope_override"]:
        print("\t", ep)
elif "protein_targets_dir" in config and len(config["protein_targets_dir"]) > 0:
    print("PROTEIN TARGETS:", )
    for f in fasta_scan(config["protein_targets_dir"]):
        print("\t", f)

fs = fasta_glob(config["proteome_search_glob"])
print("REFERENCE PROTEOME:")
print(config["proteome_search_glob"])
print("NUMBER OF FILES:")
print(len(fs))
print("FILE LIST:")
for f in fs:
    print("\t", f)
