#include "misc.h"
#include "set_params.h"
#include "io.h"
#include <ctype.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_randist.h>


void set_params_one_simul(param *params, arg *args, unsigned idx_simul){
  unsigned i, *idxs_snps_src;
  /* first we reseed the rng */
  reseed_rngs(params, args->nb_thread, "Parameters");
  /* we read all parameters from file */
  read_param_file(params, args, idx_simul);
  /* we select the SNPs we want */
  if (params->wanted_snps == NULL)
    params->wanted_snps = allocation_unsigned_vector(args->nb_snp);
  fprintf(params->simul_log, "Selecting SNPs for simulation...\n");
  if (params->total_nb_snp == args->nb_snp){
    for (i = 0; i < params->total_nb_snp; i++)
      params->wanted_snps[i] = i;
  }
  else{
    idxs_snps_src = allocation_unsigned_vector(params->total_nb_snp);
    for (i = 0; i < params->total_nb_snp; i++)
      idxs_snps_src[i] = i;
    gsl_ran_choose(params->rngs[0], params->wanted_snps, args->nb_snp, idxs_snps_src, params->total_nb_snp, sizeof(unsigned));
    free(idxs_snps_src);
  }
  /* we reseed the generators for the actual simulation */
  reseed_rngs(params, args->nb_thread, "Simulation");

}

void read_param_file(param *params, arg *args, unsigned idx_simul){
  char *filename = NULL, *buff = NULL;
  FILE * f = NULL;
  unsigned i, dummy2;
  int dummy;
  filename = allocation_char_vector(MEDIUM_BUFF_SIZE);
  buff = allocation_char_vector(MEDIUM_BUFF_SIZE);
  sprintf(filename, "./%s/simu_%u/simu_%u.par", args->prefix, idx_simul, idx_simul);
  f = safe_open(filename, "r");
  args->nb_generation = check_param_file_nb_field(f);
  /* space allocation */
  free(params->adm_Ne_simul);
  params->adm_Ne_simul = allocation_unsigned_vector(args->nb_generation+1);
  if (params->contrib_simul != NULL)
    free_double_matrix(params->contrib_simul, 2);
  params->contrib_simul = allocation_double_matrix(2, args->nb_generation+1);
  
  buff = fgets(buff, LARGE_BUFF_SIZE, f); /* skip header */
  for (i = 0; i < (args->nb_generation) + 1; i++){
    buff = fgets(buff, LARGE_BUFF_SIZE, f);
    dummy = sscanf(buff, "%u\t%u\t%lf\t%lf", &dummy2, &(params->adm_Ne_simul[i]), &(params->contrib_simul[0][i]), &(params->contrib_simul[1][i]));
    
    if (dummy != 4)
      exit_on_error(NULL, PARSING_ERROR);

    if (params->adm_Ne_simul[i] > args->max_Ne)
      exit_on_error(NULL, NE_ERROR);
    
  }
  free(filename);
  free(buff);
  fclose(f);
}

unsigned check_param_file_nb_field(FILE *f){
  char *buff = allocation_char_vector(MEDIUM_BUFF_SIZE);
  unsigned i, nb_space;
  int nb_line = 0;
  while(fgets(buff, LARGE_BUFF_SIZE, f) != NULL){
    nb_line++;
    nb_space = 0;
    for (i = 0; i < LARGE_BUFF_SIZE; i++){
      if (buff[i] == '\n')
	break;
      if (isspace(buff[i]))
	nb_space++;
    }
    if (nb_space != 3)
      exit_on_error(NULL, PARSING_ERROR);
  }
  rewind(f);
  free(buff);
  if (nb_line - 2 <= 0)
    exit_on_error(NULL, PARSING_ERROR);
  return nb_line - 2;
}

void reseed_rngs(param *params, unsigned nb_rng, char *mess){
  unsigned seed, i;
  seed = gsl_rng_get(params->rngs[0]);
  fprintf(params->simul_log, "%s seed: %u\n", mess, seed);
  gsl_rng_set(params->rngs[0], seed);
  for (i = 1; i < nb_rng; i++){
    seed = gsl_rng_get(params->rngs[0]);
    gsl_rng_set(params->rngs[i], seed);
  }
}

void reset_param(param *params, arg *args){
  fclose(params->simul_log);
  free(params->order_chroms_s1);
  free(params->order_chroms_s2);
}

void create_simul_log(arg *args, param *params, unsigned idx_simul){
  char *str = NULL;
  str = allocation_char_vector(MEDIUM_BUFF_SIZE);
  sprintf(str, "./%s/simu_%u/simu_%u.log", args->prefix, idx_simul, idx_simul);
  params->simul_log = safe_open(str, "w");
  free(str);
}


void seed_RNGs_from_file(param *params, unsigned nb_rng){
  unsigned seed, i;
  /* we get a seed */
  seed = get_seed_from_file();
  /* we initialize the first RNG with it */
  gsl_rng_set(params->rngs[0], seed);
  /* we write a new seed */
  write_next_seed(gsl_rng_get(params->rngs[0]));
  /* we initialise the remaining RNGs using the first */
  for (i = 1; i < nb_rng; i++){
    gsl_rng_set(params->rngs[i], gsl_rng_get(params->rngs[0]));
  }
}

void write_next_seed(unsigned seed){
    FILE *f = fopen("./.seed.txt", "w");
    if( f == NULL ){
        fprintf(stderr, "Cannot write next seed...\n");
        return;
    }
    fprintf(f, "%u", seed);
    fclose(f);
}


unsigned get_seed_from_file(void){
    unsigned seed;
    int read;
    FILE *f = fopen("./.seed.txt", "r");
    /* si le fichier de graine n'existe pas */
    if(f == NULL){
        fprintf(stderr, "Seeding against time...\n");
        return (unsigned) time(NULL);
    }
    read = fscanf(f, "%u", &seed);
    fclose(f);
    if(read != 1){
        fprintf(stderr, "Seeding against time...\n");
        return (unsigned) time(NULL);
    }
    return seed;
}
