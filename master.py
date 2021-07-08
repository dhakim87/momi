import json
import sqlite3
import os
from collections import defaultdict

from epitope_scanner import EpitopeScanner, EpitopeScannerCPP
from netmhciipan_wrapper import NetMHCIIpanRun
import sys

SQL_TIMEOUT = 300  # 5 minutes, rather than the default 5 seconds.  This should give time to write results from even the largets of fasta files.

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
    return sorted(paths)


def build_db():
    conn = sqlite3.connect("momi.db", timeout=SQL_TIMEOUT)
    cur = conn.cursor()
    cur.execute("CREATE TABLE IF NOT EXISTS target_epitopes (mimicked_file TEXT, mimicked_protein TEXT, mimicked_epitope TEXT, affinity TEXT)")
    cur.execute("CREATE UNIQUE INDEX IF NOT EXISTS target_epitopes_ind ON target_epitopes(mimicked_file, mimicked_protein, mimicked_epitope, affinity)")
    cur.execute("CREATE TABLE IF NOT EXISTS mimic (file TEXT, protein TEXT, offset INT, mimicked_epitope text, flanking text, epitope text, blosum real)")
    cur.execute("CREATE UNIQUE INDEX IF NOT EXISTS mimic_ind ON mimic(file,protein,offset,mimicked_epitope)")
    conn.commit()
    conn.close()


def build_target_epitopes(config):
    if "epitope_override" in config \
            and len(config["epitope_override"]) > 0 \
            and "protein_targets_dir" in config \
            and len(config["protein_targets_dir"]) > 0:
        raise Exception("Cannot specify epitope_override and protein_targets_dir simultaneously")

    if "epitope_override" in config and len(config["epitope_override"]) > 0:
        conn = sqlite3.connect("momi.db", timeout=SQL_TIMEOUT)
        cur = conn.cursor()
        results = [
            ("config", "unspecified", ep, "unspecified") for ep in config["epitope_override"]
        ]
        for r in results:
            cur.execute("INSERT OR IGNORE INTO target_epitopes (mimicked_file, mimicked_protein, mimicked_epitope, affinity) VALUES(?, ?, ?, ?)",
                        r)
        conn.commit()
        conn.close()
    else:
        for fpath in fasta_scan(config["protein_targets_dir"]):
            runner = NetMHCIIpanRun(
                config["netmhcpath"],
                config["allele"],
                config["tmpdir_prefix"],
                fpath)

            target_protein = os.path.basename(fpath)
            results = runner.run()
            conn = sqlite3.connect("momi.db", timeout=SQL_TIMEOUT)
            cur = conn.cursor()
            for protein in results:
                for core in results[protein]:
                    cur.execute("INSERT OR IGNORE INTO target_epitopes (mimicked_file, mimicked_protein, mimicked_epitope, affinity) VALUES(?, ?, ?, ?)",
                                (target_protein, protein, core, results[protein][core]))
            conn.commit()
            conn.close()


def scan_epitopes(config, job_index, num_jobs):
    conn = sqlite3.connect("momi.db", timeout=SQL_TIMEOUT)
    cur = conn.cursor()
    cur.execute("SELECT DISTINCT mimicked_epitope FROM target_epitopes")
    epitopes = [r[0] for r in cur.fetchall()]
    conn.close()

    f_index = -1
    for fpath in fasta_scan(config["proteome_search_dir"]):
        f_index += 1
        if (f_index % num_jobs) != job_index:
            continue
        print(fpath)
        if "epitope_scanner" in config and config["epitope_scanner"] != "":
            scanner_path = config['epitope_scanner']
            scanner = EpitopeScannerCPP(epitopes, fpath, config["blosum_threshold"], scanner_path)
        else:
            scanner = EpitopeScanner(epitopes, fpath, config["blosum_threshold"])
        results = scanner.run()

        conn = sqlite3.connect("momi.db", timeout=SQL_TIMEOUT)
        cur = conn.cursor()
        for result in results:
            cur.execute("INSERT OR IGNORE INTO mimic (file, protein, offset, mimicked_epitope, flanking, epitope, blosum) VALUES(?,?,?,?,?,?,?)", result)
        conn.commit()
        conn.close()


if __name__ == "__main__":
    if len(sys.argv) >= 3:
        job_index = int(sys.argv[1])
        num_jobs = int(sys.argv[2])
    else:
        job_index = 0
        num_jobs = 1

    with open("run_config.json") as config_file:
        config = json.load(config_file)
        print(config)

    build_db()
    build_target_epitopes(config)
    print("Job Index: " + str(job_index))
    print("Num Jobs: " + str(num_jobs))
    scan_epitopes(config, job_index, num_jobs)
