#include "misc.h"
#include "arguments_parser.h"
#include "io.h"
#include "set_params.h"
#include "simul.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <time.h>
#include <string.h>

/* 
 * Main function of the code
 */
int main(int argc, char **argv){
  arg *args = NULL;
  param *params = NULL;
  unsigned idx_simul;
  args = parse_arguments(argc, argv);
  params = create_params(args);
  for (idx_simul = args->idx_simul_deb; idx_simul <= args->idx_simul_fin; idx_simul++){
    fprintf(stderr, "Starting simulation %u\n", idx_simul);
    fprintf(params->general_log, "Starting simulation %u\n", idx_simul);
    create_simul_log(args, params, idx_simul);
    params->current_simul = idx_simul;
    set_params_one_simul(params, args, idx_simul);
    make_one_simul(params, args);
    reset_param(params, args);
  }
  
  free_param(params, args);
  free_arg(args);
  return EXIT_SUCCESS;
}


/* creates a param structure with general log files etc */
param *create_params(arg *args){
  param *p = NULL;
  char *name = NULL;
  name = allocation_char_vector(MEDIUM_BUFF_SIZE);
  p = allocation_param(args);
  /* then we create the log file and store it in the param structure */
  sprintf(name, "./%s/global.log", args->prefix);
  p->general_log = safe_open(name, "a");
  free(name);
  write_arguments(args, p->general_log);
  /* we allocate one RNG per thread */
  p->rngs = allocation_rng_vector(args->nb_thread);
  /* we read the input genotypes once and for all */
  read_input_genotypes(p, args);
  return p;
}

double **allocation_double_matrix(unsigned nb_row, unsigned nb_col){
  double **m = NULL;
  unsigned i;
  m = (double **) malloc(nb_row * sizeof(double*));
  check_alloc(m);
  for (i = 0; i < nb_row; i++)
    m[i] = allocation_double_vector(nb_col);
  return m;
}

double *allocation_double_vector(unsigned size){
  double *v = NULL;
  v = (double*) malloc(size * sizeof(double));
  check_alloc(v);
  return v;
}

void free_double_matrix(double **m, unsigned nb_row){
  unsigned i;
  for (i = 0; i < nb_row; i++)
    free(m[i]);
  free(m);
}

void free_unsigned_matrix(unsigned **m, unsigned nb_row){
  unsigned i;
  for (i = 0; i < nb_row; i++)
    free(m[i]);
  free(m);
}

unsigned **allocation_unsigned_matrix(unsigned nb_row, unsigned nb_col){
  unsigned **m = NULL, i;
  m = (unsigned**) malloc(nb_row * sizeof(unsigned*));
  check_alloc(m);
  for (i = 0; i < nb_row; i++)
    m[i] = allocation_unsigned_vector(nb_col);
  return m;
}

unsigned *allocation_unsigned_vector(unsigned size){
  unsigned *v = NULL;
  v = (unsigned*) malloc(size * sizeof(unsigned));
  check_alloc(v);
  return v;
}

/* Memory allocation/free and other very general functions */
void free_param(param *p, arg *args){
  unsigned i;
  fclose(p->general_log);
  for (i = 0; i < args->nb_thread; i++)
    gsl_rng_free(p->rngs[i]);
  free(p->rngs);
  free_int_matrix(p->genos_s1, p->nb_chrom_s1);
  free_int_matrix(p->genos_s2, p->nb_chrom_s2);
  free(p->adm_Ne_simul);
  free_double_matrix(p->contrib_simul, 2);
  free(p->wanted_snps);
  free(p->microsat_mut_rates);
  free(p->SNI_mut_rates);
  free(p->geom_params);
  free(p);
}

gsl_rng **allocation_rng_vector(unsigned size){
  unsigned i;
  gsl_rng **rngs = NULL;
  const gsl_rng_type *T = gsl_rng_mt19937;
  rngs = (gsl_rng**) malloc(size * sizeof(gsl_rng*));
  check_alloc(rngs);
  for (i = 0; i < size; i++){
    rngs[i] = gsl_rng_alloc(T);
    check_alloc(rngs[i]);
  }
  /* we seed the first RNG against time */
  gsl_rng_set(rngs[0], (unsigned) time(NULL));
  return rngs;
}

/* opens and file and captures errors */
FILE *safe_open(char *name, char *mode){
  FILE *f = NULL;
  f = fopen(name, mode);
  if (f == NULL)
    exit_on_error(name, FOPEN_ERROR);
  return f;
}

/* allocation of param structure */
param *allocation_param(arg *args){
  param *p = NULL;
  p = (param*) malloc(sizeof(param));
  check_alloc(p);
  p->nb_chrom_s1 = 0;
  p->nb_chrom_s2 = 0;
  p->current_adm_Ne = 0;
  p->prev_adm_Ne = 0;
  p->prev_adm_Ne = 0;
  p->genos_s1 = NULL;
  p->genos_s2 = NULL;
  p->genos_adm_old = NULL;
  p->genos_adm_new = NULL;
  p->wanted_snps = NULL;
  p->order_chroms_s1 = NULL;
  p->order_chroms_s2 = NULL;
  p->order_chroms_adm = NULL;
  p->adm_Ne_simul = NULL;
  p->contrib_simul = NULL;
  p->microsat_mut_rates = allocation_double_vector(args->nb_snp);
  p->SNI_mut_rates = allocation_double_vector(args->nb_snp);
  p->geom_params = allocation_double_vector(args->nb_snp);
  return p;
}

couple **allocation_couple_vector(unsigned size){
  couple **v = NULL;
  unsigned i;
  v = (couple**) malloc(size * sizeof(couple*));
  check_alloc(v);
  for (i = 0; i < size; i++){
    v[i] = (couple*) malloc(size * sizeof(couple));
    check_alloc(v[i]);
    v[i]->genos1 = allocation_int_ptr_vector(2);
    v[i]->genos2 = allocation_int_ptr_vector(2);
    v[i]->chroms_names = allocation_char_matrix(4, SMALL_BUFF_SIZE);
  }
  return v;
}

void free_couple_vector(couple **v, unsigned size){
  unsigned i;
  for (i = 0; i < size; i++){
    free_char_matrix(v[i]->chroms_names, 4);
    free(v[i]->genos1);
    free(v[i]->genos2);
    free(v[i]);
  }
  free(v);
  return;
}

char **allocation_char_ptr_vector(unsigned size){
  char **v = NULL;
  unsigned i;
  v = (char **) malloc(size * sizeof(char*));
  check_alloc(v);
  for (i = 0; i < size; i++)
    v[i] = NULL;
  return v;
}

int **allocation_int_ptr_vector(unsigned size){
  int **v = NULL;
  unsigned i;
  v = (int **) malloc(size * sizeof(int*));
  check_alloc(v);
  for (i = 0; i < size; i++)
    v[i] = NULL;
  return v;
}

couple_param *allocation_couple_param(void){
  couple_param *p = NULL;
  p = (couple_param*) malloc(sizeof(couple_param));
  check_alloc(p);
  p->name_pop1 = allocation_char_vector(SMALL_BUFF_SIZE);
  p->name_pop2 = allocation_char_vector(SMALL_BUFF_SIZE);
  reset_couple_param(p);
  return p;
}

void free_couple_param(couple_param *p){
  free(p->name_pop1);
  free(p->name_pop2);
  free(p);
  return;
}

/* allocation of arg structure */
arg *allocation_arg(void){
  arg *a = NULL;
  a = (arg*) malloc(sizeof(arg));
  check_alloc(a);
  a->help = 0;
  a->nb_snp = 50000;
  a->save_data = 0;
  a->nb_generation = 5;
  a->idx_simul_deb = 1;
  a->idx_simul_fin = 1;
  a->nb_thread = 1;
  a->sample_size_s1 = 100;
  a->sample_size_s2 = 100;
  a->sample_size_adm = 100;
  a->input_path = allocation_char_vector(MEDIUM_BUFF_SIZE);
  a->prefix = allocation_char_vector(MEDIUM_BUFF_SIZE);
  a->input_path[0] = '\0';
  a->prefix[0] = '\0';
  return a;
}

pthread_t *allocation_pthread_t_vector(unsigned size){
  pthread_t *v = NULL;
  v = (pthread_t*) malloc(size * sizeof(pthread_t));
  check_alloc(v);
  return v;
}

thread_datast *allocation_thread_datast_vector(unsigned size){
  thread_datast *v = NULL;
  v = (thread_datast*) malloc(size * sizeof(thread_datast));
  check_alloc(v);
  return v;
}

void free_arg(arg *a){
  free(a->input_path);
  free(a->prefix);
  free(a);
}

char **allocation_char_matrix(unsigned nb_line, unsigned nb_col){
  char **m = NULL;
  unsigned i;
  m = (char**) malloc(nb_line * sizeof(char*));
  check_alloc(m);
  for (i = 0; i < nb_line; i++)
    m[i] = allocation_char_vector(nb_col);
  return m;
}

void free_char_matrix(char **m, unsigned nb_line){
  unsigned i;
  for (i = 0; i < nb_line; i++)
    free(m[i]);
  free(m);
}

char *allocation_char_vector(unsigned size){
  char *v = NULL;
  v = (char*) malloc(size * sizeof(char));
  check_alloc(v);
  return v;
}

int **allocation_int_matrix(unsigned nb_line, unsigned nb_col){
  int **m = NULL;
  unsigned i;
  m = (int**) malloc(nb_line * sizeof(int*));
  check_alloc(m);
  for (i = 0; i < nb_line; i++)
    m[i] = allocation_int_vector(nb_col);
  return m;
}

void free_int_matrix(int **m, unsigned nb_line){
  unsigned i;
  for (i = 0; i < nb_line; i++)
    free(m[i]);
  free(m);
}

int *allocation_int_vector(unsigned size){
  int *v = NULL;
  v = (int*) malloc(size * sizeof(int));
  check_alloc(v);
  return v;
}


/* checks a pointer is not NULL
 * calls 'exit_on_error' if it is the case
 */
void check_alloc(void *p){
  if (p == NULL)
    exit_on_error(NULL, MEMALLOC_ERROR);
  return;
}

/* prints usage message and exists the program */
void print_usage(char *progname){
  fprintf(stderr, "Usage: %s OPTIONS\n", progname);
  fprintf(stderr, "\t--help\n");
  fprintf(stderr, "\t\tDisplay this help and exit\n");
  fprintf(stderr, "\t--save-data\n");
  fprintf(stderr, "\t\tSave sampled genotypes at last generation (default=False)\n");
  fprintf(stderr, "\t--sampling S1/ADM/S2\n");
  fprintf(stderr, "\t\tNumber of individuals to sample for summary statistics computations (default=100/100/100)\n");
  fprintf(stderr, "\t--nb-snp N\n");
  fprintf(stderr, "\t\tNumber of SNP to simulate (default=50000)\n");
  fprintf(stderr, "\t--prefix\n");
  fprintf(stderr, "\t\tThe path where we will find all simulations directories\n");
  fprintf(stderr, "\t--nb-simul N\n");
  fprintf(stderr, "\t\tNumber of simulations to perform (default=1)\n");
  fprintf(stderr, "\t--nb-thread N\n");
  fprintf(stderr, "\t\tNumber of parallel processes to run (default=1)\n");
  fprintf(stderr, "\t--input-path\n");
  fprintf(stderr, "\t\tPath to the input file containing genotypes for source populations 1 and 2 in Arlequin format (fsc output)\n");
  exit_on_error(NULL, NO_CODE);
}

/* Prints an error message and exists the program */
void exit_on_error(char *arg, unsigned error_code){
  switch (error_code){
  case NO_CODE :
    break;
  case MEMALLOC_ERROR :
    fprintf(stderr, "Memory allocation error. Exiting\n");
    break;
  case FOPEN_ERROR :
    fprintf(stderr, "Cannot open file %s. Exiting\n", arg);
    break;
  case SSIZE_ERROR :
    fprintf(stderr, "%s", arg);
    break;
  case PARSING_ERROR :
    fprintf(stderr, "An error occured while parsing input file. Exiting\n");
    break;
  case NB_SNP_ERROR :
    fprintf(stderr, "Not enough SNP in input file to simulate the number of SNPs required. Exiting\n");
    break;
  case BIG_ROUNDING_ERROR :
    fprintf(stderr, "A large rounding error occured, simulation is compromised. Exiting\n");
    break;
  case PTHREAD_ERROR :
    fprintf(stderr, "An error occured with the pthread library. Exiting\n");
    break;
  }
  exit(EXIT_FAILURE);
}

