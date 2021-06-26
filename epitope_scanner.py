import os
import sys
from uuid import uuid4


def readBlosum62():
    dict = {}
    with open("blosum62X.txt") as f:
        data = f.read().splitlines()
        aas = data[0].split()

        i = 0
        for row in data[1:]:
            vals = row.split()[1:]
            for i in range(len(vals)):
                score = int(vals[i])
                colAA = aas[i]
                rowAA = row.split()[0]
                dict[(rowAA, colAA)] = score
        return dict

BLOSUM62 = readBlosum62()


class EpitopeScanner:
    def __init__(self, epitopes, fpath, blosum_thresh):
        self.epitopes = epitopes
        self.fpath = fpath
        self.fname = os.path.basename(self.fpath)

        self.thresh = blosum_thresh
        self.flanking = 6
        self.epitope_len = len(self.epitopes[0])

    def _calc_blosum(self, putative_mimic):
        scores = {}
        for core_epitope in self.epitopes:
            sum = 0
            if len(core_epitope) != len(putative_mimic):
                raise Exception("Epitope length mismatch")
            for i in range(len(core_epitope)):
                key = (core_epitope[i], putative_mimic[i])
                if key in BLOSUM62:
                    sum += BLOSUM62[key]
                else:
                    print(core_epitope)
                    print(putative_mimic)
                    raise Exception("Unknown Amino Acid Match: ", key)
            if sum >= self.thresh:
                scores[core_epitope] = sum
        return scores

    def _scan(self, protein_name, protein_sequence, out_arr):
        protein_sequence = ("^"*self.flanking) + protein_sequence + ("$"*self.flanking)
        for window_start in range(self.flanking, len(protein_sequence) - self.flanking - self.epitope_len):
            window_end = window_start + self.epitope_len
            window = protein_sequence[window_start: window_end]

            blosum_hits = self._calc_blosum(window)
            # For every match
            # - new uuid
            # - file_name
            # - protein_name
            # - flanking window
            # - mimic
            # - mimicked_epitope
            # - blosum score
            if len(blosum_hits) > 0:
                flanked_window = protein_sequence[window_start-self.flanking: window_end+self.flanking]
            for mimicked_epitope in blosum_hits:
                score = blosum_hits[mimicked_epitope]
                out_arr.append(
                    [str(uuid4()), self.fname, protein_name, flanked_window, window, mimicked_epitope, score]
                )

    def run(self):
        out_arr = []
        active = []
        protein_name = ""
        i = 0
        with open(self.fpath) as fasta:
            for line in fasta:
                if line.startswith(">"):
                    i += 1
                    print(i)
                    self._scan(protein_name, "".join(active), out_arr)
                    protein_name = line[1:-1]
                    active = []
                else:
                    active.append(line[:-1])
            if protein_name != "":
                self._scan(protein_name, "".join(active), out_arr)

        return out_arr