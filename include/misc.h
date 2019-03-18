
#ifndef MISC_H
#define MISC_H

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <pthread.h>

/* definition of general codes */
#define NO_CODE 0
#define MEMALLOC_ERROR 1
#define CREATE_DIR_ERROR 2
#define FOPEN_ERROR 3
#define SSIZE_ERROR 4
#define PARSING_ERROR 5
#define NB_SNP_ERROR 6
#define BIG_ROUNDING_ERROR 7
#define PTHREAD_ERROR 8

#define MOD_CON 0
#define MOD_INC 1
#define MOD_DEC 2
#define MOD_ALL 3
#define PLOT_NONE 0
#define PLOT_LAST 1
#define PLOT_ALL 2
#define TREND_NONE 0
#define TREND_EXP 1

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
  unsigned adm_Ne0;
  unsigned adm_final_Ne_min;
  unsigned adm_final_Ne_max;
  unsigned *adm_Ne_simul;
  unsigned *order_chroms_s1;
  unsigned *order_chroms_s2;
  unsigned *order_chroms_adm;
  double contrib_0_min[2];
  double contrib_0_max[2];
  double contrib_1_min[2];
  double contrib_1_max[2];
  double contrib_final_min[2];
  double contrib_final_max[2];
  double **contrib_simul;
  unsigned *wanted_snps;
  char **genos_adm_old;
  char **genos_adm_new;
}param;


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

/* /\* structure containing infos for an admixed individual *\/ */
/* typedef struct indiv{ */
/*   unsigned *parents; */
/*   char **chroms; */
/* }indiv; */


/* structure containing arguments from the command line */
typedef struct arg{
  int help;
  int save_data;
  int force_rewrite;
  unsigned draw_plots;
  int draw_preview;
  unsigned nb_snp;
  unsigned nb_generation;
  unsigned idx_simul_deb;
  unsigned idx_simul_fin;
  unsigned nb_thread;
  unsigned sample_size_s1;
  unsigned sample_size_s2;
  unsigned sample_size_adm;
  unsigned adm_Ne0;
  unsigned adm_Ne_model;
  unsigned adm_final_Ne_default;
  unsigned adm_final_Ne_min;
  unsigned adm_final_Ne_max;
  unsigned contrib_0_default[2];
  double contrib_0[2];
  unsigned contrib_model[2];
  unsigned contrib_trend[2];
  double contrib_trend_min[2];
  double contrib_trend_max[2];
  unsigned contrib_1_default[2];
  double contrib_1_min[2];
  double contrib_1_max[2];
  unsigned contrib_final_default[2];
  double contrib_final_min[2];
  double contrib_final_max[2];
  char *input_path;
  char *prefix;
}arg;


/* general functions definitions */
void exit_on_error(char *arg, unsigned error_code);
void check_alloc(void *p);
arg *allocation_arg(void);
void free_arg(arg *a);
char *allocation_char_vector(unsigned size);
void print_usage(char *progname);
param *allocation_param(arg *args);
param *create_params(arg *args);
void create_directory(char *name, unsigned force_rewrite);
FILE *safe_open(char *name, char *mode);
void free_param(param *p, arg *args);
gsl_rng **allocation_rng_vector(unsigned size);
char **allocation_char_matrix(unsigned nb_line, unsigned nb_col);
void free_char_matrix(char **m, unsigned nb_line);
unsigned *allocation_unsigned_vector(unsigned size);
double **allocation_double_matrix(unsigned nb_row, unsigned nb_col);
void free_double_matrix(double **m, unsigned nb_row);
double *allocation_double_vector(unsigned size);
/* indiv **allocation_indiv_vector(unsigned size, unsigned nb_snp); */
/* indiv *allocation_indiv(unsigned nb_snp); */
/* void free_indiv_vector(indiv **v, unsigned size); */
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