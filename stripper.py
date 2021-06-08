import sys

# NetMHCIIPan chokes on peptides shorter than its search window.  We remove
# fasta proteins with length < 15.  

def write_protein(out_file, protein_lines):
    # We check against 16 because lines from file have newlines
    if len(protein_lines) < 2 or len(protein_lines[1]) < 16:
        return
    for pl in protein_lines:
        out_file.write(pl)


protein_lines = []
with open(sys.argv[1] + ".stripped", "w") as out_file:
    with open(sys.argv[1]) as fasta:
        for line in fasta:
            if line.startswith(">"):
                write_protein(out_file, protein_lines)
                protein_lines = [line]
            else:
                protein_lines.append(line)
            
        write_protein(out_file, protein_lines)

