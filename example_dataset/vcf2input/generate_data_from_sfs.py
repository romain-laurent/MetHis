#!/usr/bin/env python

import sys
import numpy as np

def get_sfs(path1, path2) :
    f1 = open(path1)
    f2 = open(path2)
    f1.next()                   # headers
    f2.next()
    sfs = []
    for line1 in f1 :
        line2 = f2.next()
        line1 = line1.split()
        line2 = line2.split()
        if line1[0] != line2[0] or line1[1] != line2[1] :
            print >> sys.stderr, 'Erreur de position'
            sys.exit(1)
        line1 = line1[-2:]
        line2 = line2[-2:]
        line1 = [i.split(':') for i in line1]
        line2 = [i.split(':') for i in line2]
        if line1[0][0] == line2[0][0] :
            freqs = [float(line1[0][-1]), float(line2[0][-1])]
            sfs.append(freqs)
            continue
        if line1[0][0] == line2[1][0] :
            freqs = [float(line1[0][-1]), float(line2[1][-1])]
            sfs.append(freqs)
            continue
        print >> sys.stderr, 'Erreur d\'alleles'
        sys.exit()
    f1.close()
    f2.close()
    return sfs

def generate_data(sfs, idx, wanted_nb_chrom) :
    data = []
    for i in xrange(len(sfs)) :
        f = sfs[i][idx]
        chroms = np.random.choice(['0','1'], wanted_nb_chrom, p=[f,1-f])
        data.append(''.join(list(chroms)))
    return data

def write_data(data, idx_pop, nb_chrom, path) :
    f = open(path, 'a')
    print >> f, '{'
    for i in xrange(nb_chrom) :
        new_line = ['{0}_{1}'.format(idx_pop, i+1), ' ', '1', ' ']
        for j in xrange(len(data)) :
            new_line.append(data[j][i])
        print >> f, ''.join(new_line)
    print >> f, '}'
    f.close()

def main() :
    path_sfs1 = 'sfs_s1.frq'
    path_sfs2 = 'sfs_s2.frq'
    path_out = 'MetHis_input.arp'
    f = open(path_out, 'w')
    f.close()
    wanted_nb_chrom = 20000
    sfs_both_pops = get_sfs(path_sfs1, path_sfs2)
    data = generate_data(sfs_both_pops, 0, wanted_nb_chrom)
    write_data(data, 1, wanted_nb_chrom, path_out)
    data = generate_data(sfs_both_pops, 1, wanted_nb_chrom)
    write_data(data, 2, wanted_nb_chrom, path_out)

if __name__ == '__main__' :
    main()
