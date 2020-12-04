#ifndef SUMSTATS_H
#define SUMSTATS_H

#define NB_DIM_MDS 2
#define NB_POINTS_KDE 512	/* default R value, this number MUST be a power of 2 */
#define NB_PERCENTILE 11

double compute_a_Fst(double n_bar, double n_C, double ssquare, double p_bar, double h_bar);
double compute_b_Fst(double n_bar, double ssquare, double p_bar, double h_bar);
double compute_ssquare_Fst(unsigned n1, unsigned n2, double p1, double p2, double n_bar, double p_bar);
double compute_bar_stat_Fst(unsigned ssize1, unsigned ssize2, double tilde1, double tilde2, double n_bar);
void initialize_thread_datast_freqs(thread_datast_freqs *data, arg *args);
void finalize_inbreeding_computations(thread_datast_freqs *threads_data, sumstats_datast *sumstats, arg *args);
void *compute_allelic_frequencies_snp_subset(void *ptr);
void free_freqs_datast(freqs_datast *f);
freqs_datast *compute_allelic_frequencies_local(char **genotypes, arg *args, sumstats_datast *sumstats);
freqs_datast *allocation_freqs_datast(arg *args);
void initialize_sumstats(sumstats_datast *s);
double *compute_kords(double hi, double lo, double bw);
double *compute_bin_dist(double *data, double lo, double hi, unsigned size);
double **compute_kernel_density_estimate(double *data, unsigned size);
void compute_ASD_stats(double **asd, arg *args, sumstats_datast *sumstats);
double nrd0(double *data, unsigned size);
double gaussian_kernel(double x);
double compute_distribution_mode(double *data, unsigned size);
double compute_dot_product(double *v1, double *v2);
double *compute_vector_two_points(double *source, double *dest);
double compute_dist_two_points(double *p1, double *p2);
double *compute_centroid(double **coords, unsigned idx_min, unsigned idx_max);
double *compute_admixture_proportions(double **mds, arg *args, unsigned nb_indiv);
double *compute_admixture_angles(double **mds, arg *args, unsigned nb_indiv);
double **compute_mds(double **dists, unsigned nb_indiv);
void *compute_asd_snp_subset(void *ptr);
double **compute_asd_matrix_local(param *params, arg *args, char **genotypes);
thread_datast_asd *allocation_thread_datast_asd_vector(unsigned nb_thread);
void compute_all_sumstats(param *params, arg *args, unsigned **wanted_samples);
char **simplify_genotypes(param *params, arg *args, unsigned **wanted_samples);
void *fill_genotype_matrix(void *ptr);
thread_datast_fill_genotypes *allocation_thread_datast_fill_genotypes_vector(unsigned nb_thread);
char recode_geno(char geno1, char geno2);

#endif
