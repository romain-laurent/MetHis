#!/usr/bin/env python2

import sys
import argparse
import random
import os.path
import os

global max_nb_try
max_nb_try = 10000

# creates parser for command line arguments and parses them
def parse_args() :
    parser = argparse.ArgumentParser(description='Generates parameters for \'methis\' simulations.')
    parser.add_argument('-S', '--nb_simulation',
                        type=int,
                        default=1,
                        help='Number of parameters sets to create',
                        dest='nb_simulation')
    parser.add_argument('-N', '--nb_generation',
                        type=int,
                        required=True,
                        help='Number of generations to simulate',
                        dest='nb_generation')
    parser.add_argument('-P', '--prefix',
                        required=True,
                        help='Prefix for parameters files names',
                        dest='prefix')
    parser.add_argument('--Ne',
                        required=True,
                        dest='Ne',
                        help='Description of admixed population Ne evolution over time. See documentation for help about syntax')
    parser.add_argument('--contrib_s1',
                        required=True,
                        dest='contrib_s1',
                        help='Description of s1 population contribution to admixed population over time. See documentation for help about syntax')
    parser.add_argument('--contrib_s2',
                        required=True,
                        dest='contrib_s2',
                        help='Description of s2 population contribution to admixed population over time. See documentation for help about syntax')
    parser.add_argument('--force-rewrite',
                        default=False,
                        action='store_true',
                        dest='force_rewrite',
                        help='Should we rewrite a file if it already exists? Default=False')
    
                        
    values = parser.parse_args()
    return values.nb_simulation, values.nb_generation, values.prefix, values.Ne, values.contrib_s1, values.contrib_s2, values.force_rewrite

# prints error message and exits if parsing of Ne fails
def exit_on_Ne_str_parsing_fail(Ne_str) :
    print >> sys.stderr, 'Failed parsing Ne parameter: {0}'.format(Ne_str)
    print >> sys.stderr, 'Exiting'
    sys.exit(1)

# prints the interpretation of Ne argument
def print_Ne_interpretation(Ne0, scheme, Ne_range) :
    print >> sys.stderr, 'Ne argument interpreted as:'
    print >> sys.stderr, '\tNe at founding generation: {0}'.format(Ne0)
    print >> sys.stderr, '\tNe variation scheme: {0}'.format(scheme)
    print >> sys.stderr, '\tNe final range: {0}'.format(Ne_range)
    
# parse the Ne string from command line
def parse_Ne_str(Ne_str) :
    Ne_splitted = Ne_str.split('/')
    if len(Ne_splitted) != 3 :
        exit_on_Ne_str_parsing_fail(Ne_str)
    try :
        Ne0 = int(Ne_splitted[0])
    except :
        exit_on_Ne_str_parsing_fail(Ne_str)
    if Ne_splitted[1] not in ('Con','Inc','Dec','All') :
        exit_on_Ne_str_parsing_fail(Ne_str)
    if Ne_splitted[2] == 'default' :
        Ne_range = (10, 5000)
    else :
        try :
            tmp = Ne_splitted[2].split('-')
            Ne_range = int(tmp[0]), int(tmp[1])
            Ne_range = tuple(sorted(Ne_range))
        except :
            exit_on_Ne_str_parsing_fail(Ne_str)
    print_Ne_interpretation(Ne0, Ne_splitted[1], Ne_range)
    return Ne0, Ne_splitted[1], Ne_range

# exits in case of contrib string parsing fail
def exit_on_contrib_parsing_fail(contrib_str, contrib_id) :
    print >> sys.stderr, 'Failed parsing contrib_s{0} parameter: {1}'.format(contrib_id, contrib_str)
    print >> sys.stderr, 'Exiting'
    sys.exit(1)

def print_contrib_interpretation(contrib, contrib_id) :
    print >> sys.stderr, 'Contribution {0} interpreted as:'.format(contrib_id)
    print >> sys.stderr, '\tFounding contribution: {0}'.format(contrib[0])
    print >> sys.stderr, '\tContribution scheme: {0}'.format(contrib[1])
    print >> sys.stderr, '\tInitial contribution range: {0}'.format(contrib[2])
    print >> sys.stderr, '\tFinal contribution range: {0}'.format(contrib[3])

# parses a 'Pulse' contribution scheme
def parse_pulse_contrib(contrib_str, contrib_id) :
    contrib_splitted = contrib_str.split('/')
    if contrib_splitted[0] == 'default' :
        c0 = 'default'
    else :
        try :
            c0 = float(contrib_splitted[0])
        except :
            exit_on_contrib_parsing_fail(contrib_str, contrib_id)
    try :
        nb_pulse = int(contrib_splitted[2])
        if len(contrib_splitted) != nb_pulse + 3 :
            raise ValueError
    except :
        exit_on_contrib_parsing_fail(contrib_str, contrib_id)
    contrib = [c0, 'Pulse', nb_pulse]
    for i in xrange(nb_pulse) :
        tmp = contrib_splitted[i+3]
        if tmp == 'default' :
            contrib.append((0.0,1.0))
        else :
            try :
                tmp = tmp.split('-')
                if len(tmp) != 2 :
                    raise ValueError
                tmp = tuple(sorted([float(tmp[0]), float(tmp[1])]))
                contrib.append(tmp)
            except :
                exit_on_contrib_parsing_fail(contrib_str, contrib_id)
    print >> sys.stderr, 'Contribution {0} interpreted as:'.format(contrib_id)
    print >> sys.stderr, '\tFounding contribution: {0}'.format(contrib[0])
    print >> sys.stderr, '\tContribution scheme: Pulse'
    for i in xrange(contrib[2]) :
        print >> sys.stderr, '\tPulse {0} range: {1}'.format(i+1, contrib[3+i])
    return contrib
   
# parses a contribution string from the command line
def parse_contrib_str(contrib_str, contrib_id) :
    # we make a special case for the 'Pulse' scheme
    if 'Pulse' in contrib_str :
        contrib = parse_pulse_contrib(contrib_str, contrib_id)
        return contrib
    # for other schemes
    contrib_splitted = contrib_str.split('/')
    if len(contrib_splitted) != 4 :
        exit_on_contrib_parsing_fail(contrib_str, contrib_id)
    if contrib_splitted[0] == 'default' :
        contrib_0 = 'default'
    else :
        try :
            contrib_0 = float(contrib_splitted[0])
        except :
            exit_on_contrib_parsing_fail(contrib_str, contrib_id)
    if contrib_splitted[1] not in ('Con','Inc','Dec','All', 'Pulse') :
        exit_on_contrib_parsing_fail(contrib_str, contrib_id)
    if contrib_splitted[2] == 'default' :
        contrib_range_init = (0.0, 1.0)
    else :
        try :
            tmp = contrib_splitted[2].split('-')
            if len(tmp) != 2 :
                raise ValueError
            contrib_range_init = float(tmp[0]), float(tmp[1])
            contrib_range_init = tuple(sorted(contrib_range_init))
        except :
            exit_on_contrib_parsing_fail(contrib_str, contrib_id)
    if contrib_splitted[3] == 'default' :
        contrib_range_final = (0.0, 1.0)
    else :
        try :
            tmp = contrib_splitted[3].split('-')
            if len(tmp) != 2 :
                raise ValueError
            contrib_range_final = float(tmp[0]), float(tmp[1])
            contrib_range_final = tuple(sorted(contrib_range_final))
        except :
            exit_on_contrib_parsing_fail(contrib_str, contrib_id)
    contrib = (contrib_0, contrib_splitted[1], contrib_range_init, contrib_range_final)
    print_contrib_interpretation(contrib, contrib_id)
    return contrib

def exit_on_Ne_inconsistency() :
    print >> sys.stderr, 'ERROR: Ne pattern is inconsistent. Exiting'
    sys.exit(1)

# checks that Ne pattern makes sense
def make_Ne_sanity_check(Ne) :
    Ne0 = Ne[0]
    scheme = Ne[1]
    min_final, max_final = Ne[2]
    if scheme == 'Inc' and max_final < Ne0 :
        exit_on_Ne_inconsistency()
    if scheme == 'Dec' and min_final > Ne0 :
        exit_on_Ne_inconsistency()
    return Ne0, scheme, (min_final, max_final)

def exit_on_contrib_inconsistency(contrib_id) :
    print >> sys.stderr, 'ERROR: Contribution {0} pattern is inconsistent. Exiting'.format(contrib_id)
    sys.exit(1)
        
# checks that contribution patterns make sense
def make_contrib_sanity_check(contrib, contrib_id) :
    c0 = contrib[0]
    scheme = contrib[1]
    if scheme == 'Pulse' :
        return contrib
    c_init_min, c_init_max = contrib[2]
    c_final_min, c_final_max = contrib[3]
    if scheme == 'Inc' and c_init_min > c_final_max :
        exit_on_contrib_inconsistency(contrib_id)
    if scheme == 'Inc' :
        c_init_max = min(c_init_max, c_final_max)
        c_final_min = max(c_final_min, c_init_min)
        print >> sys.stderr, 'WARNING: Adjusting ranges for contribution {0}'.format(contrib_id)
        print >> sys.stderr, '\tInitial contribution range: {0}-{1}'.format(c_init_min, c_init_max)
        print >> sys.stderr, '\tFinal contribution range: {0}-{1}'.format(c_final_min, c_final_max)
    if scheme == 'Dec' and c_init_max < c_final_min :
        exit_on_contrib_inconsistency(contrib_id)
    if scheme == 'Dec' :
        c_init_min = max(c_init_min, c_final_min)
        c_final_max = min(c_final_max, c_init_max)
        print >> sys.stderr, 'WARNING: Adjusting ranges for contribution {0}'.format(contrib_id)
        print >> sys.stderr, '\tInitial contribution range: {0}-{1}'.format(c_init_min, c_init_max)
        print >> sys.stderr, '\tFinal contribution range: {0}-{1}'.format(c_final_min, c_final_max)
    if scheme == 'Con' :
        if c_init_max < c_final_min or c_init_min > c_final_max :
            exit_on_contrib_inconsistency(contrib_id)
        c_init_min = max(c_init_min, c_final_min)
        c_init_max = min(c_init_max, c_final_max)
        c_final_min = c_init_min
        c_final_max = c_init_max
        print >> sys.stderr, 'WARNING: Adjusting ranges for contribution {0}'.format(contrib_id)
        print >> sys.stderr, '\tInitial contribution range: {0}-{1}'.format(c_init_min, c_init_max)
        print >> sys.stderr, '\tFinal contribution range: {0}-{1}'.format(c_final_min, c_final_max)
    if c_init_min < 0 or c_init_max > 1 or c_final_min < 0 or c_final_max > 1 :
        exit_on_contrib_inconsistency(contrib_id)
    if not isinstance(c0, str) :
        if c0 < 0 or c0 > 1 :
            exit_on_contrib_inconsistency(contrib_id)
    return c0, scheme, (c_init_min, c_init_max), (c_final_min, c_final_max)

# exits if contributions patterns are not compatible
def exit_on_compatibility_failure() :
    print >> sys.stderr, 'ERROR: contributions patterns are not compatible. Exiting'
    sys.exit(1)

# checks that contributions patterns are compatible
def check_contributions_consistency(p1, p2) :
    # initial contributions should sum to 1
    if p1[0] != 'default' and p2[0] != 'default' :
        if p1[0] + p2[0] != 1 :
            exit_on_compatibility_failure()
    if p1[1] == 'Pulse' or p2[1] == 'Pulse' :
        return
    if p1[2][0] + p2[2][0] > 1 :
        exit_on_compatibility_failure()
    if p1[3][0] + p2[3][0] > 1 :
        exit_on_compatibility_failure()

# generates a growing or decreasing pattern using Functional response type 2 equation (https://en.wikipedia.org/wiki/Functional_response)
# the math implementation creates decreasing values, which we reverse if we need an increasing pattern
def generate_non_constant_pattern(x1, x2, y1, y2, increase) :
    if x2 == x1 :
        print >> sys.stderr, 'ERROR: cannot generate non constant pattern over 1 generation.\nExiting'
        sys.exit(1)
    u = random.uniform(0, 0.5)
    a = u*u / (1 - 2*u)
    values = []
    for i in xrange(x1, x2 + 1) :
        tmp_i = float(i)
        new_val = (y1-y2) * a * (1-((tmp_i-x1)/(x2-x1))) / (a+((tmp_i-x1)/(x2-x1))) + y2
        values.append(new_val)
    if increase :
        values.reverse()
    return u, values
        
# Generates Ne vector from Ne pattern
def generate_Nes(pattern, nb_generation) :
    Ne0 = pattern[0]
    scheme = pattern[1]
    final_min, final_max = pattern[2]
    new_Nes = []
    # we use the 'failed' procedure in case some of the schemes of 'All' are not possible with the given ranges
    failed = True
    tmp_final = 0
    while failed :
        if scheme == 'All' :
            tmp_scheme = random.choice(['Con','Inc','Dec'])
        else :
            tmp_scheme = scheme
        if tmp_scheme == 'Con' :
            tmp = random.randint(final_min, final_max)
            new_Nes = [Ne0]
            new_Nes.extend([tmp] * nb_generation)
            u = 'NA'
            failed = False
        if tmp_scheme == 'Inc' :
            tmp_final = random.randint(final_min, final_max)
            if tmp_final < Ne0 :
                continue
            u, Nes = generate_non_constant_pattern(1, nb_generation+1, tmp_final, Ne0, True)
            new_Nes.extend(Nes)
            new_Nes = [int(round(j)) for j in new_Nes]
            failed = False
        if tmp_scheme == 'Dec' :
            tmp_final = random.randint(final_min, final_max)
            if tmp_final > Ne0 :
                continue
            u, Nes = generate_non_constant_pattern(1, nb_generation+1, Ne0, tmp_final, False)
            new_Nes.extend(Nes)
            new_Nes = [int(round(j)) for j in new_Nes]
            failed = False
    return u, new_Nes

# Generates a pulse contribution
def generate_pulse_contribution(pattern, nb_generation) :
    c0 = pattern[0]
    nb_pulse = pattern[2]
    generations = range(1, nb_generation + 1)
    pulses_times = sorted(random.sample(generations, nb_pulse))
    new_contribs = [c0]
    new_contribs.extend([0] * nb_generation)
    intensities = []
    for i in xrange(nb_pulse) :
        i_min, i_max = pattern[i+3]
        intensities.append(random.uniform(i_min, i_max))
        new_contribs[pulses_times[i]] = intensities[-1]
    return (pulses_times, intensities), new_contribs

# Generates contribution vector from contribution pattern
def generate_contribution(pattern, nb_generation) :
    c0 = pattern[0]
    scheme = pattern[1]
    # we make a special case for Pulse scheme
    if scheme == 'Pulse' :
        u, new_contribs = generate_pulse_contribution(pattern, nb_generation)
        return u, new_contribs
    init_min, init_max = pattern[2]
    final_min, final_max = pattern[3]
    new_contribs = [c0]
    # we use the 'failed' procedure in case some of the schemes of 'All' are not possible with the given ranges
    failed = True
    while failed :
        if scheme == 'All' :
            tmp_scheme = random.choice(['Con','Inc','Dec'])
        else :
            tmp_scheme = scheme
        if tmp_scheme == 'Con' :
            tmp = random.uniform(final_min, final_max)
            new_contribs.extend([tmp] * nb_generation)
            u = 'NA'
            failed = False
        if tmp_scheme == 'Inc' :
            tmp_final = 1.0
            tmp_init = 1.0
            while tmp_final < 3*tmp_init :
                tmp_final = random.uniform(final_min, final_max)
                tmp_init = random.uniform(init_min, init_max)
            if tmp_final < tmp_init :
                continue
            u, contrib = generate_non_constant_pattern(1, nb_generation, tmp_final, tmp_init, True)
            new_contribs.extend(contrib)
            failed = False
        if tmp_scheme == 'Dec' :
            tmp_final = 1.0
            tmp_init = 1.0
            while tmp_final*3 > tmp_init :
                tmp_init = random.uniform(init_min, init_max)
                tmp_final = random.uniform(final_min, final_max)
            if tmp_final > tmp_init :
                continue
            u, contrib = generate_non_constant_pattern(1, nb_generation, tmp_init, tmp_final, False)
            new_contribs.extend(contrib)
            failed = False
    return u, new_contribs

# exits if we fail to genererate contributions too many times
def exit_on_nb_try(max_nb_try) :
    print >> sys.stderr, 'Failed to generate valid contributions pattern after {0} tries. Exiting'.format(max_nb_try)
    sys.exit(1)
    
# Once new contributions are available, checks that at no point them sum to more than 1
def check_contribs_compatibility(c1, c2) :
    for i in xrange(len(c1)) :
        if c1[i] + c2[i] > 1 :
            return False
    return True

# exits if the number of pulses is bigger than the number of generations to simulate
def exit_on_nb_pulse_error(contrib_id) :
    print >> sys.stderr, 'ERROR: too much pulses for contribution {0}.\nExiting'.format(contrib_id)
    sys.exit(1)

# special function for generating c0s in case user used the 'default' value
def generate_c0s(contribs_s1, contribs_s2) :
    if isinstance(contribs_s1[0], float) and isinstance(contribs_s2[0], float) :
        return contribs_s1, contribs_s2
    if isinstance(contribs_s1[0], str) and isinstance(contribs_s2[0], str) :
        contribs_s1[0] = random.random()
        contribs_s2[0] = 1 - contribs_s1[0]
        return contribs_s1, contribs_s2
    if isinstance(contribs_s1[0], str) :
        contribs_s1[0] = 1 - contribs_s2[0]
    else :
        contribs_s2[0] = 1 - contribs_s1[0]
    return contribs_s1, contribs_s2

# Prints real parameters we want to estimate latter
def print_real_parameters(prefix, idx_simu, u1, contribs_s1, u2, contribs_s2, uNe, Nes) :
    header = ['Ne.0', 'Ne.1', 'Ne.N', 'Ne.u']
    params = [str(Nes[0]), str(Nes[1]), str(Nes[-1]), str(uNe)]
    header.append('s1.0')
    params.append(str(contribs_s1[0]))
    if isinstance(u1, tuple) :
        times, intensities = u1
        nb_pulse = len(times)
        for i in xrange(nb_pulse) :
            header.append('t1.{0}'.format(i+1))
            params.append(str(times[i]))
            header.append('i1.{0}'.format(i+1))
            params.append(str(intensities[i]))
    else :
        header.extend(['s1.1', 's1.N', 's1.u'])
        params.extend([str(contribs_s1[1]), str(contribs_s1[-1]), str(u1)])
    
    header.append('s2.0')
    params.append(str(contribs_s2[0]))
    if isinstance(u2, tuple) :
        times, intensities = u2
        nb_pulse = len(times)
        for i in xrange(nb_pulse) :
            header.append('t2.{0}'.format(i+1))
            params.append(str(times[i]))
            header.append('i2.{0}'.format(i+1))
            params.append(str(intensities[i]))
    else :
        header.extend(['s2.1', 's2.N', 's2.u'])
        params.extend([str(contribs_s2[1]), str(contribs_s2[-1]), str(u2)])
    
    f = open('{0}/simu_{1}/simu_{1}.txt'.format(prefix, idx_simu), 'w')
    print >> f, '\t'.join(header)
    print >> f, '\t'.join(params)
    f.close()



# main function
def main() :
    global max_nb_try
    # parsing command line
    nb_simulation, nb_generation, prefix, Ne_str, contrib_s1_str, contrib_s2_str, force_rewrite = parse_args()

    # parsing patterns
    Ne_pattern = parse_Ne_str(Ne_str)
    s1_pattern = parse_contrib_str(contrib_s1_str, 1)
    s2_pattern = parse_contrib_str(contrib_s2_str, 2)
    
    # sanity checks + adjusting range
    Ne_pattern = make_Ne_sanity_check(Ne_pattern)
    s1_pattern = make_contrib_sanity_check(s1_pattern, 1)
    s2_pattern = make_contrib_sanity_check(s2_pattern, 2)
    
    check_contributions_consistency(s1_pattern, s2_pattern)
    if s1_pattern[1] == 'Pulse' and len(s1_pattern) > nb_generation + 3 :
        exit_on_nb_pulse_error(1)
    if s2_pattern[1] == 'Pulse' and len(s2_pattern) > nb_generation + 3 :
        exit_on_nb_pulse_error(2)
    try :
        os.mkdir('{0}'.format(prefix))
    except OSError :
        if not force_rewrite :
            print >> sys.stderr, 'Directory {0} already exists. Exiting'.format(prefix)
            sys.exit(1)
    # starting parameters generation
    for i in xrange(1, nb_simulation+1) :
        uNe, new_Nes = generate_Nes(Ne_pattern, nb_generation)
        contribs_compatible = False
        nb_try = 0
        while not contribs_compatible :
            if nb_try > max_nb_try :
                exit_on_nb_try(max_nb_try)
            nb_try += 1
            u1, contribs_s1 = generate_contribution(s1_pattern, nb_generation)
            u2, contribs_s2 = generate_contribution(s2_pattern, nb_generation)
            contribs_s1, contribs_s2 = generate_c0s(contribs_s1, contribs_s2)
            contribs_compatible = check_contribs_compatibility(contribs_s1, contribs_s2)

        # printing result to file
        filename = '{0}/simu_{1}/simu_{1}.par'.format(prefix, i)
        if os.path.isfile(filename) and force_rewrite == False :
            print >> sys.stderr, 'File {0} exists. Doing nothing. Check the force-rewrite argument'.format(filename)
            continue
        try :
            os.mkdir('{0}/simu_{1}'.format(prefix, i))
        except :
            pass
        f = open(filename, 'w')
        print >> f, 'g\tNe\tc1\tc2'
        for j in xrange(nb_generation+1) :
            print >> f, '\t'.join([str(j), str(new_Nes[j]), str(contribs_s1[j]), str(contribs_s2[j])])
        f.close()
        # printing real parameters we want to estimate latter
        print_real_parameters(prefix, i, u1, contribs_s1, u2, contribs_s2, uNe, new_Nes)

if __name__ == '__main__' :
    main()
