
#ifndef MISC_H
#define MISC_H

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <pthread.h>

/* definition of general codes */
#define NO_CODE 0
#define MEMALLOC_ERROR 1
#define FOPEN_ERROR 3
#define SSIZE_ERROR 4
#define PARSING_ERROR 5
#define NB_SNP_ERROR 6
#define BIG_ROUNDING_ERROR 7
#define PTHREAD_ERROR 8
#define NE_ERROR 9


#define FLAG_ADM 0
#define FLAG_S1 1
#define FLAG_S2 2


/* definition of general constants */
#define SMALL_BUFF_SIZE 64
#define MEDIUM_BUFF_SIZE 1024
#define LARGE_BUFF_SIZE 1048576
#define MAX_NE 10000
#define MIN_NE 20

/* definition of general macros */
#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) < (b) ? (b) : (a))

#define get_idx(x) ((x) == '1' ? 0 : 1)
#define get_id(x) ((x) == 0 ? '1' : '2')


#define PROG_NAME "MetHis"


/* structure containing general parameters + parameters for one simulation */
typedef struct param{
  FILE *general_log;
  FILE *simul_log;
  gsl_rng **rngs;
  unsigned nb_chrom_s1;
  unsigned nb_chrom_s2;
  char **genos_s1;
  char **genos_s2;
  unsigned total_nb_snp;
  unsigned current_simul;
  unsigned current_adm_Ne;
  unsigned prev_adm_Ne;
  unsigned *adm_Ne_simul;
  unsigned *order_chroms_s1;
  unsigned *order_chroms_s2;
  unsigned *order_chroms_adm;
  double **contrib_simul;
  unsigned *wanted_snps;
  char **genos_adm_old;
  char **genos_adm_new;
}param;


/* structure containing arguments from the command line */
typedef struct arg{
  int help;
  int save_data;
  unsigned nb_snp;
  unsigned max_Ne;
  unsigned nb_generation;
  unsigned idx_simul_deb;
  unsigned idx_simul_fin;
  unsigned nb_thread;
  unsigned sample_size_s1;
  unsigned sample_size_s2;
  unsigned sample_size_adm;
  char *input_path;
  char *prefix;
}arg;


/* structure used to store thread information for genotype conversion */
typedef struct thread_datast_fill_genotypes{
  unsigned idx_min;
  unsigned idx_max;
  param *params;
  arg *args;
  unsigned **wanted_samples;
  char **genotypes;
} thread_datast_fill_genotypes;

/* structure used to store thread information when computing allele sharing distances */
typedef struct thread_datast_asd{
  unsigned idx_min;
  unsigned idx_max;
  unsigned nb_indiv;
  char **genotypes;
  double **tmp_asd;
} thread_datast_asd;

/* structure used to store allele frequencies */
typedef struct freqs_datast{
  double *freqs_adm;
  double *freqs_s1;
  double *freqs_s2;
  double *hets_adm;
  double *hets_s1;
  double *hets_s2;
} freqs_datast;

/* structure used to store thread information for computation of statistics based on allele frequencies */
typedef struct thread_datast_freqs{
  char **genotypes;
  freqs_datast *freqs;
  arg *args;
  unsigned idx_deb;
  unsigned idx_fin;
  double accu_num_f3;
  double accu_denom_f3;
  unsigned *nb_het_sites_adm;
  unsigned *nb_het_sites_s1;
  unsigned *nb_het_sites_s2;
  unsigned nb_sites_considered_F_adm;
  unsigned nb_sites_considered_F_s1;
  unsigned nb_sites_considered_F_s2;
  double expected_hom_adm;
  double expected_hom_s1;
  double expected_hom_s2;
  double n_bar_adm_s1;
  double n_bar_adm_s2;
  double n_bar_s1_s2;
  double n_C_adm_s1;
  double n_C_adm_s2;
  double n_C_s1_s2;
  double accu_num_fst_adm_s1;
  double accu_denom_fst_adm_s1;
  double accu_num_fst_adm_s2;
  double accu_denom_fst_adm_s2;
  double accu_num_fst_s1_s2;
  double accu_denom_fst_s1_s2;
}thread_datast_freqs;

/* structure used to store summary statistics */
typedef struct sumstats_datast{
  double adm_prop_mean;
  double adm_prop_sd;
  double adm_prop_skew;
  double adm_prop_kurt;
  double adm_prop_mode;
  double *percentiles;
  double adm_ang_mean;
  double adm_ang_sd;
  double adm_ang_skew;
  double adm_ang_kurt;
  double adm_ang_mode;
  double *ang_percentiles;
  
  double mean_asd_adm;
  double var_asd_adm;
  double mean_asd_s1;
  double var_asd_s1;
  double mean_asd_s2;
  double var_asd_s2;
  double mean_asd_adm_s1;
  double var_asd_adm_s1;
  double mean_asd_adm_s2;
  double var_asd_adm_s2;
  double mean_asd_s1_s2;
  double var_asd_s1_s2;
  double f3;
  double mean_het_adm;
  double var_het_adm;
  double mean_het_s1;
  double var_het_s1;
  double mean_het_s2;
  double var_het_s2;
  double mean_F_adm;
  double var_F_adm;
  double mean_F_s1;
  double var_F_s1;
  double mean_F_s2;
  double var_F_s2;
  double fst_adm_s1;
  double fst_adm_s2;
  double fst_s1_s2;
}sumstats_datast;



/* structure containing informations about one given couple */
typedef struct couple{
  char **genos1;
  char **genos2;
  char **chroms_names;
}couple;

/* structure contraining infos needed to run one thread */
typedef struct thread_datast{
  unsigned idx_min;
  unsigned idx_max;
  unsigned nb_snp;
  couple **couples;
  param *params;
  gsl_rng *rng;
}thread_datast;


/* structure containing information for couples creations */
typedef struct couple_param{
  char **genos_pop1;
  char **genos_pop2;
  unsigned nb_indiv_pop1;
  unsigned nb_indiv_pop2;
  char *name_pop1;
  char *name_pop2;
  unsigned *order_pop1;
  unsigned *order_pop2;
}couple_param;

/* general functions definitions */
void exit_on_error(char *arg, unsigned error_code);
void check_alloc(void *p);
arg *allocation_arg(void);
void free_arg(arg *a);
char *allocation_char_vector(unsigned size);
void print_usage(char *progname);
param *allocation_param(arg *args);
param *create_params(arg *args);
FILE *safe_open(char *name, char *mode);
void free_param(param *p, arg *args);
gsl_rng **allocation_rng_vector(unsigned size);
char **allocation_char_matrix(unsigned nb_line, unsigned nb_col);
void free_char_matrix(char **m, unsigned nb_line);
unsigned *allocation_unsigned_vector(unsigned size);
double **allocation_double_matrix(unsigned nb_row, unsigned nb_col);
void free_double_matrix(double **m, unsigned nb_row);
double *allocation_double_vector(unsigned size);
couple_param *allocation_couple_param(void);
void free_couple_param(couple_param *p);
couple **allocation_couple_vector(unsigned size);
char **allocation_char_ptr_vector(unsigned size);
void free_couple_vector(couple **v, unsigned size);
thread_datast *allocation_thread_datast_vector(unsigned size);
pthread_t *allocation_pthread_t_vector(unsigned size);
unsigned **allocation_unsigned_matrix(unsigned nb_row, unsigned nb_col);
void free_unsigned_matrix(unsigned **m, unsigned nb_row);

#endif
