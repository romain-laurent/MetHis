
#ifndef SIMUL_H
#define SIMUL_H

#define S1S1 0
#define S1H 1
#define HH 2
#define S2H 3
#define S2S2 4
#define S1S2 5

void make_one_simul(param *params, arg *args);
void count_nb_couple_each_class(unsigned *nb_couples, unsigned Ne, double p1, double p2, gsl_rng *rng);
unsigned *initialize_chrom_order(unsigned nb_chrom);
void reset_couple_param(couple_param *p);
void set_couple_param(couple_param *p, unsigned couple_type, param *params);
void create_one_couple(couple *to_create, couple_param *couple_params, gsl_rng *r);
void make_multithreaded_reproduction(param *params, couple **couples, unsigned nb_thread, unsigned nb_snp);
void *make_reproduction(void *ptr);
unsigned **sample_individuals(param *params, arg *args, couple **couples);
unsigned adm_indivs_share_chromosomes(couple *couple1, couple *couple2);
unsigned adm_shares_chromosome_with_source(char *chrom1, char *chrom2, unsigned **wanted, couple **couples, unsigned nb_adm);
void compute_asd_matrix(arg *args, param *params, unsigned idx_gen);
void remove_temporary_files(arg *args, param *params, unsigned idx_gen);
void remove_one_file(char *name);
void compute_allelic_frequencies(param *params, arg *args, unsigned idx_gen);
void compute_fst(param *params, arg *args, unsigned idx_gen);
void compute_inbreeding(param *params, arg *args, unsigned idx_gen);
void finalize_stats_computation(param *params, arg *args, unsigned idx_gen);

#endif
