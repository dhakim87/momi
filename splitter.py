import sys 

num_parts = int(sys.argv[2])

fds = []
for i in range(num_parts):
    fds.append(open(sys.argv[1] + "." + str(i), 'w'))

j = -1
out_file = None
with open(sys.argv[1]) as fasta:
    for line in fasta:
        if line.startswith(">"):
            j += 1
            j = j % num_parts
            out_file = fds[j]
        out_file.write(line)

for fd in fds:
    fd.close()

