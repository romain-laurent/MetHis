#!/usr/bin/env python3

import random

nb_wanted = 60000

snps = set()
f = open('indep_SNPs.prune.in')
for line in f :
    snps.add(line.strip())
f.close()

snps = list(snps)
snps = random.sample(snps, nb_wanted)

for i in snps :
    print(i)
