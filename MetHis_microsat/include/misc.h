
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
  int **genos_s1;
  int **genos_s2;
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
  int **genos_adm_old;
  int **genos_adm_new;
  double *microsat_mut_rates;
  double *SNI_mut_rates;
  double *geom_params;
}param;


/* structure containing informations about one given couple */
typedef struct couple{
  int **genos1;
  int **genos2;
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
  int **genos_pop1;
  int **genos_pop2;
  unsigned nb_indiv_pop1;
  unsigned nb_indiv_pop2;
  char *name_pop1;
  char *name_pop2;
  unsigned *order_pop1;
  unsigned *order_pop2;
}couple_param;

/* structure containing arguments from the command line */
typedef struct arg{
  int help;
  int save_data;
  unsigned nb_snp;
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

int **allocation_int_matrix(unsigned nb_row, unsigned nb_col);
void free_int_matrix(int **m, unsigned nb_row);
int *allocation_int_vector(unsigned size);
int **allocation_int_ptr_vector(unsigned size);

#endif
