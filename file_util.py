import os
import glob
import gzip
import lzma


# Transparently handles gz and xz files
def open_fasta(fpath, mode):
    if fpath.endswith(".fasta.gz"):
        return gzip.open(fpath, mode=mode, encoding='utf-8')
    elif fpath.endswith(".fasta.xz"):
        return lzma.open(fpath, mode=mode, encoding='utf-8')
    else:
        return open(fpath, mode=mode)


def _check_file(fpath):
    if not os.path.isfile(fpath):
        return False
    tup = os.path.splitext(fpath)
    if len(tup) < 2:
        return False
    target_protein, extension = tup
    if extension not in [".fasta", ".fasta.gz", ".fsa", ".faa", ".fna", ".ffn", ".faa", ".frn", ".fa"] \
            and not fpath.endswith(".fasta.gz")\
            and not fpath.endswith(".fasta.xz"):  # Ugh bioinformatics.
        return False
    return True


def fasta_scan(dir):
    paths = []
    for target_protein_file in os.listdir(dir):
        fpath = os.path.join(dir, target_protein_file)
        if _check_file(fpath):
            paths.append(fpath)
    return sorted(paths)


def fasta_glob(glob_str):
    paths = []
    for fpath in glob.glob(glob_str):
        if _check_file(fpath):
            paths.append(fpath)
    return sorted(paths)
