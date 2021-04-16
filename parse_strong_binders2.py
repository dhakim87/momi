import sys

cores = set([])
for line in sys.stdin:
    ss = line.split()
    # Example split line: ['322', 'DRB1_1501', 'TGVTRLYHAARLRYK', '3', 'TRLYHAARL', '0.993', 'EGS31750.1', '0.568872', '1.35', 'NA', '<=SB']
    
    pos, allele, seq, of, core, core_rel, identity, score_el, rank_el, exp_bind, bind_level = ss

    #print(core, identity, bind_level)
    cores.add((core, identity))

for core in sorted(list(cores)):
    print(",".join(core))
