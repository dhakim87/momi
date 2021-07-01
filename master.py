import json
import sqlite3
import os
from collections import defaultdict

from epitope_scanner import EpitopeScanner, EpitopeScannerCPP
from netmhciipan_wrapper import NetMHCIIpanRun


def fasta_scan(dir):
    paths = []
    for target_protein_file in os.listdir(dir):
        fname = target_protein_file
        fpath = os.path.join(dir, target_protein_file)
        if not os.path.isfile(fpath):
            continue
        tup = os.path.splitext(fname)
        if len(tup) < 2:
            continue
        target_protein, extension = tup
        if extension not in [".fasta", ".fsa", ".faa", ".fna", ".ffn", ".faa", ".frn", ".fa"]:  #Ugh bioinformatics.
            continue
        paths.append(fpath)
    return paths


def build_db():
    conn = sqlite3.connect("momi.db")
    cur = conn.cursor()
    cur.execute("CREATE TABLE IF NOT EXISTS target_epitopes (mimicked_file TEXT, mimicked_protein TEXT, mimicked_epitope TEXT, affinity TEXT)")
    cur.execute("CREATE UNIQUE INDEX IF NOT EXISTS target_epitopes_ind ON target_epitopes(mimicked_file, mimicked_protein, mimicked_epitope, affinity)")
    cur.execute("CREATE TABLE IF NOT EXISTS mimic (file TEXT, protein TEXT, offset INT, mimicked_epitope text, flanking text, epitope text, blosum real)")
    cur.execute("CREATE UNIQUE INDEX IF NOT EXISTS mimic_ind ON mimic(file,protein,offset,mimicked_epitope)")
    conn.commit()
    conn.close()


def build_target_epitopes(config):
    for fpath in fasta_scan(config["protein_targets_dir"]):
        runner = NetMHCIIpanRun(
            config["netmhcpath"],
            config["allele"],
            config["tmpdir_prefix"],
            fpath)

        target_protein = os.path.basename(fpath)
        results = runner.run()
        conn = sqlite3.connect("momi.db")
        cur = conn.cursor()
        for protein in results:
            for core in results[protein]:
                cur.execute("INSERT OR IGNORE INTO target_epitopes (mimicked_file, mimicked_protein, mimicked_epitope, affinity) VALUES(?, ?, ?, ?)",
                            (target_protein, protein, core, results[protein][core]))
        conn.commit()
        conn.close()


def scan_epitopes(config):
    conn = sqlite3.connect("momi.db")
    cur = conn.cursor()
    cur.execute("SELECT DISTINCT mimicked_epitope FROM target_epitopes")
    epitopes = [r[0] for r in cur.fetchall()]
    conn.close()
    for fpath in fasta_scan(config["proteome_search_dir"]):
        # scanner = EpitopeScanner(epitopes, fpath, config["blosum_threshold"])
        scanner = EpitopeScannerCPP(epitopes, fpath, config["blosum_threshold"], "./a.out")
        results = scanner.run()

        conn = sqlite3.connect("momi.db")
        cur = conn.cursor()
        for result in results:
            cur.execute("INSERT OR IGNORE INTO mimic (file, protein, offset, mimicked_epitope, flanking, epitope, blosum) VALUES(?,?,?,?,?,?,?)", result)
        conn.commit()
        conn.close()


if __name__ == "__main__":
    with open("run_config.json") as config_file:
        config = json.load(config_file)
        print(config)

    build_db()
    build_target_epitopes(config)
    scan_epitopes(config)