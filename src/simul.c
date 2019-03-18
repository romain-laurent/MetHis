#include "misc.h"
#include "simul.h"
#include "io.h"
#include <string.h>
#include <math.h>
#include <gsl/gsl_randist.h>


void make_one_simul(param *params, arg *args){
  unsigned idx_gen, i, j, nb_done;
  unsigned nb_couples[6];
  couple_param *couple_params = NULL;
  couple **couples = NULL;
  unsigned **wanted_samples = NULL;
  /* we initialise the order_chroms vectors */
  params->order_chroms_s1 = initialize_chrom_order(params->nb_chrom_s1);
  params->order_chroms_s2 = initialize_chrom_order(params->nb_chrom_s2);
  couple_params = allocation_couple_param();
  for (idx_gen = 0; idx_gen <= args->nb_generation; idx_gen++){
    params->current_adm_Ne = params->adm_Ne_simul[idx_gen];
    if (idx_gen != 0)
      params->prev_adm_Ne = params->adm_Ne_simul[idx_gen-1];
    fprintf(stderr, "\tStarting generation %u...\n", idx_gen);
    fprintf(params->simul_log, "\tStarting generation %u...\n", idx_gen);
    /* we count the number of couples of each class contributing to adm pop */
    count_nb_couple_each_class(nb_couples, params->current_adm_Ne, params->contrib_simul[0][idx_gen], params->contrib_simul[1][idx_gen], params->rngs[0]);
    /* we shuffle the order of chromosomes in the source populations */
    gsl_ran_shuffle(params->rngs[0], params->order_chroms_s1, params->nb_chrom_s1, sizeof(unsigned));
    gsl_ran_shuffle(params->rngs[0], params->order_chroms_s2, params->nb_chrom_s2, sizeof(unsigned));
    /* now we proceed with reproduction */
    params->genos_adm_new = allocation_char_matrix(2*params->current_adm_Ne, args->nb_snp);
    if (idx_gen != 0){
      params->order_chroms_adm = initialize_chrom_order(2*params->prev_adm_Ne);
    }
    /* REPRODUCTION */
    /* first we create couples */
    nb_done = 0;
    couples = allocation_couple_vector(params->current_adm_Ne);
    for (i = 0; i < 6; i++){
      if (nb_couples[i] == 0)
	continue;
      set_couple_param(couple_params, i, params);
      for (j = 0; j < nb_couples[i]; j++){
	create_one_couple(couples[nb_done], couple_params, params->rngs[0]);
	nb_done++;
      }
      reset_couple_param(couple_params);
    }
    /* now we shuffle the couples (this is important for sampling different kind of individuals later) */
    gsl_ran_shuffle(params->rngs[0], couples, params->current_adm_Ne, sizeof(couple*));
    /* and now we make individuals reproduce */
    make_multithreaded_reproduction(params, couples, args->nb_thread, args->nb_snp);
    /* we sample individuals if we need them */
    if (args->draw_plots == PLOT_ALL || idx_gen == args->nb_generation){
      /* sample individuals */
      wanted_samples = sample_individuals(params, args, couples);

      /* now we compute summary statistics */
      /* write genotypes to tped file */
      write_tped_file(params, args, wanted_samples, idx_gen);
      /* we compute ASD matrix using external software */
            compute_asd_matrix(args, params, idx_gen);
      /* if we want to draw plot */
      if (args->draw_plots != PLOT_NONE)
	draw_MDS_plots(args, params, idx_gen);
      /* write genotypes to vcf file */
      write_vcf_file(params, args, wanted_samples, idx_gen);
      compute_allelic_frequencies(params, args, idx_gen);

      /* /\* create input for ADMIXTOOLS *\/ */
      /* create_admixtools_inputs(params, args, idx_gen); */
      /* /\* /\\* compute f3 statistic *\\/ *\/ */
      /* compute_f3_stat(params, args); */

      compute_fst(params, args, idx_gen);
      compute_inbreeding(params, args, idx_gen);
      finalize_stats_computation(params, args, idx_gen);
      /* we remove temporary files */
      remove_temporary_files(args, params, idx_gen);
      /* at the last generation, we compute summary statistics anyway */
      free_unsigned_matrix(wanted_samples, args->sample_size_s1 + args->sample_size_s2 + args->sample_size_adm);
    }
    /* a little bit of cleanup */
    if (idx_gen != 0){
      free_char_matrix(params->genos_adm_old, 2*params->prev_adm_Ne);
      free(params->order_chroms_adm);
    }
    params->genos_adm_old = params->genos_adm_new;
    free_couple_vector(couples, params->current_adm_Ne);
  }
  free_char_matrix(params->genos_adm_old, 2*params->current_adm_Ne);
  free_couple_param(couple_params);
  return;
}

void finalize_stats_computation(param *params, arg *args, unsigned idx_gen){
  char *command = NULL;
  int dummy;
  fprintf(stderr, "Finalizing summary statistics computations...\n");
  fprintf(params->simul_log, "Finalizing summary statistics computations...\n");
  command = allocation_char_vector(LARGE_BUFF_SIZE);
  sprintf(command, "/usr/bin/Rscript R_codes/finalize_sumstats.R %s %u %u", args->prefix, params->current_simul, idx_gen);
  dummy = system(command);
  if (dummy != 0){
    fprintf(stderr, "It seems something went wrong during final sumstats computation\n");
    fprintf(params->simul_log, "It seems something went wrong during final sumstats computation\n");
  }
  free(command);
}

void compute_inbreeding(param *params, arg *args, unsigned idx_gen){
  char *command = NULL, *tmp = NULL, *tmp2 = NULL;
  int dummy;
  fprintf(stderr, "Computing inbreeding coefficients...\n");
  fprintf(params->simul_log, "Computing inbreeding coefficients...\n");
  command = allocation_char_vector(LARGE_BUFF_SIZE);
  tmp = allocation_char_vector(MEDIUM_BUFF_SIZE);
  tmp2 = allocation_char_vector(MEDIUM_BUFF_SIZE);
  sprintf(tmp, "%s/simu_%u/simu_%u_g%u", args->prefix, params->current_simul, params->current_simul, idx_gen);
  sprintf(tmp2, "%s/simu_%u", args->prefix, params->current_simul);
  sprintf(command, "%s --vcf %s.vcf --keep %s/indivs_adm.txt --het --out %s/adm 2> /dev/null", VCFTOOLS_PATH, tmp, tmp2, tmp2);
  dummy = system(command);
  if (dummy != 0){
    fprintf(stderr, "It seems something went wrong during inbreeding coefficients computation\n");
    fprintf(params->simul_log, "It seems something went wrong during inbreeding coefficients computation\n");
  }
  sprintf(command, "%s --vcf %s.vcf --keep %s/indivs_s1.txt --het --out %s/s1 2> /dev/null", VCFTOOLS_PATH, tmp, tmp2, tmp2);
  dummy = system(command);
  if (dummy != 0){
    fprintf(stderr, "It seems something went wrong during inbreeding coefficients computation\n");
    fprintf(params->simul_log, "It seems something went wrong during inbreeding coefficients computation\n");
  } 
  sprintf(command, "%s --vcf %s.vcf --keep %s/indivs_s2.txt --het --out %s/s2 2> /dev/null", VCFTOOLS_PATH, tmp, tmp2, tmp2);
  dummy = system(command);
  if (dummy != 0){
    fprintf(stderr, "It seems something went wrong during inbreeding coefficients computation\n");
    fprintf(params->simul_log, "It seems something went wrong during inbreeding coefficients computation\n");
  }
  
  free(command);
  free(tmp);
  free(tmp2);
}


void compute_fst(param *params, arg *args, unsigned idx_gen){
  char *command = NULL, *tmp = NULL, *tmp2 = NULL;
  int dummy;
  fprintf(stderr, "Computing Fst...\n");
  fprintf(params->simul_log, "Computing Fst...\n");
  command = allocation_char_vector(LARGE_BUFF_SIZE);
  tmp = allocation_char_vector(MEDIUM_BUFF_SIZE);
  tmp2 = allocation_char_vector(MEDIUM_BUFF_SIZE);
  sprintf(tmp, "%s/simu_%u/simu_%u_g%u", args->prefix, params->current_simul, params->current_simul, idx_gen);
  sprintf(tmp2, "%s/simu_%u", args->prefix, params->current_simul);
  sprintf(command, "%s --vcf %s.vcf --weir-fst-pop %s/indivs_adm.txt --weir-fst-pop %s/indivs_s1.txt --out %s 2>&1 | grep weighted | awk '{print $7}' > %s/result.fst", VCFTOOLS_PATH, tmp, tmp2, tmp2, tmp, tmp2);
  dummy = system(command);
  if (dummy != 0){
    fprintf(stderr, "It seems something went wrong during Fst computation\n");
    fprintf(params->simul_log, "It seems something went wrong during Fst computation\n");
  }
  sprintf(command, "%s --vcf %s.vcf --weir-fst-pop %s/indivs_adm.txt --weir-fst-pop %s/indivs_s2.txt --out %s 2>&1 | grep weighted | awk '{print $7}' >> %s/result.fst", VCFTOOLS_PATH, tmp, tmp2, tmp2, tmp, tmp2);
  dummy = system(command);
  if (dummy != 0){
    fprintf(stderr, "It seems something went wrong during Fst computation\n");
    fprintf(params->simul_log, "It seems something went wrong during Fst computation\n");
  }
  sprintf(command, "%s --vcf %s.vcf --weir-fst-pop %s/indivs_s1.txt --weir-fst-pop %s/indivs_s2.txt --out %s 2>&1 | grep weighted | awk '{print $7}' >> %s/result.fst", VCFTOOLS_PATH, tmp, tmp2, tmp2, tmp, tmp2);
  dummy = system(command);
  if (dummy != 0){
    fprintf(stderr, "It seems something went wrong during Fst computation\n");
    fprintf(params->simul_log, "It seems something went wrong during Fst computation\n");
  }
  free(command);
  free(tmp);
  free(tmp2);
}

void compute_allelic_frequencies(param *params, arg *args, unsigned idx_gen){
  char *command = NULL, *tmp = NULL, *tmp2 = NULL;
  int dummy;
  fprintf(stderr, "Computing allelic frequencies...\n");
  fprintf(params->simul_log, "Computing allelic frequencies...\n");
  command = allocation_char_vector(LARGE_BUFF_SIZE);
  tmp = allocation_char_vector(MEDIUM_BUFF_SIZE);
  tmp2 = allocation_char_vector(MEDIUM_BUFF_SIZE);
  sprintf(tmp, "%s/simu_%u/simu_%u_g%u", args->prefix, params->current_simul, params->current_simul, idx_gen);
  sprintf(tmp2, "%s/simu_%u", args->prefix, params->current_simul);
  sprintf(command, "%s --vcf %s.vcf --keep %s/indivs_adm.txt --freq2 --out %s_adm 2> /dev/null", VCFTOOLS_PATH, tmp, tmp2, tmp);
  dummy = system(command);
  if (dummy != 0){
    fprintf(stderr, "It seems something went wrong during allelic frequencies computation\n");
    fprintf(params->simul_log, "It seems something went wrong during allelic frequencies computation\n");
  }
  sprintf(command, "%s --vcf %s.vcf --keep %s/indivs_s1.txt --freq2 --out %s_s1 2> /dev/null", VCFTOOLS_PATH, tmp, tmp2, tmp);
  dummy = system(command);
  if (dummy != 0){
    fprintf(stderr, "It seems something went wrong during allelic frequencies computation\n");
    fprintf(params->simul_log, "It seems something went wrong during allelic frequencies computation\n");
  }
  sprintf(command, "%s --vcf %s.vcf --keep %s/indivs_s2.txt --freq2 --out %s_s2 2> /dev/null", VCFTOOLS_PATH, tmp, tmp2, tmp);
  dummy = system(command);
  if (dummy != 0){
    fprintf(stderr, "It seems something went wrong during allelic frequencies computation\n");
    fprintf(params->simul_log, "It seems something went wrong during allelic frequencies computation\n");
  }
  free(command);
  free(tmp);
  free(tmp2);
}

/* removes files used to compute statistics / draw MDS */
void remove_temporary_files(arg *args, param *params, unsigned idx_gen){
  char *tmp = NULL, *tmp2 = NULL;
  int dummy;
  fprintf(stderr, "Removing temporary files...\n");
  fprintf(params->simul_log, "Removing temporary files...\n");
  tmp = allocation_char_vector(MEDIUM_BUFF_SIZE);
  tmp2 = allocation_char_vector(MEDIUM_BUFF_SIZE);
  sprintf(tmp, "%s/simu_%u/simu_%u_g%u", args->prefix, params->current_simul, params->current_simul, idx_gen);
  sprintf(tmp2, "%s.asd.dist", tmp);
  remove_one_file(tmp2);
  sprintf(tmp2, "%s.error", tmp);
  remove_one_file(tmp2);
  sprintf(tmp2, "%s.log", tmp);
  remove_one_file(tmp2);
  sprintf(tmp2, "%s_adm.frq", tmp);
  remove_one_file(tmp2);
  sprintf(tmp2, "%s_s1.frq", tmp);
  remove_one_file(tmp2);
  sprintf(tmp2, "%s_s2.frq", tmp);
  remove_one_file(tmp2);
  sprintf(tmp2, "%s.weir.fst", tmp);
  remove_one_file(tmp2);
  sprintf(tmp2, "%s.vcf", tmp);
  if (args->save_data){
    sprintf(tmp, "bgzip -f %s", tmp2);
    dummy = system(tmp);
    if (dummy != 0){
      fprintf(stderr, "It seems something went wrong during data compression\n");
      fprintf(params->simul_log, "It seems something went wrong during data compression\n");
    }
  }
  else {
    remove_one_file(tmp2);
  }
  sprintf(tmp, "%s/simu_%u", args->prefix, params->current_simul);
  sprintf(tmp2, "%s/indivs_adm.txt", tmp);
  remove_one_file(tmp2);
  sprintf(tmp2, "%s/indivs_s1.txt", tmp);
  remove_one_file(tmp2);
  sprintf(tmp2, "%s/indivs_s2.txt", tmp);
  remove_one_file(tmp2);
  sprintf(tmp2, "%s/adm.het", tmp);
  remove_one_file(tmp2);
  sprintf(tmp2, "%s/s1.het", tmp);
  remove_one_file(tmp2);
  sprintf(tmp2, "%s/s2.het", tmp);
  remove_one_file(tmp2);
  sprintf(tmp2, "%s/result.fst", tmp);
  remove_one_file(tmp2);
  
  /* if we do not want to save the data */
  sprintf(tmp, "%s/simu_%u/simu_%u_g%u", args->prefix, params->current_simul, params->current_simul, idx_gen);
  sprintf(tmp2, "%s.tped", tmp);
  remove_one_file(tmp2);
  sprintf(tmp2, "%s.tfam", tmp);
  remove_one_file(tmp2);
  free(tmp);
  free(tmp2);
}

void remove_one_file(char *name){
  int dummy;
  dummy = remove(name);
  if (dummy != 0){
    fprintf(stderr, "An error occured while removing \"%s\"\n", name);
  }
}

void draw_MDS_plots(arg *args, param *params, unsigned idx_gen){
  char *command = NULL, *tmp = NULL;
  int dummy;
  fprintf(stderr, "Drawing MDS plots...\n");
  fprintf(params->simul_log, "Drawing MDS plots...\n");
  command = allocation_char_vector(LARGE_BUFF_SIZE);
  tmp = allocation_char_vector(MEDIUM_BUFF_SIZE);
  sprintf(tmp, "%s/simu_%u/simu_%u_g%u", args->prefix, params->current_simul, params->current_simul, idx_gen);
  sprintf(command, "Rscript R_codes/plot_MDS.R %s > /dev/null 2> /dev/null", tmp);
  dummy = system(command);
  if (dummy != 0){
    fprintf(stderr, "It seems something went wrong during MDS plotting\n");
    fprintf(params->simul_log, "It seems something went wrong during MDS plotting\n");
  }
  free(command);
  free(tmp);
}

void compute_asd_matrix(arg *args, param *params, unsigned idx_gen){
  char *command = NULL, *tmp = NULL;
  int dummy;
  fprintf(stderr, "Computing ASD matrix...\n");
  fprintf(params->simul_log, "Computing ASD matrix...\n");
  command = allocation_char_vector(LARGE_BUFF_SIZE);
  tmp = allocation_char_vector(MEDIUM_BUFF_SIZE);
  sprintf(tmp, "%s/simu_%u/simu_%u_g%u", args->prefix, params->current_simul, params->current_simul, idx_gen);
  sprintf(command, "%s --tfam %s.tfam --tped %s.tped --biallelic --out %s --threads %u > /dev/null 2> /dev/null", ASD_PATH, tmp, tmp, tmp, args->nb_thread);
  dummy = system(command);
  if (dummy != 0){
    fprintf(stderr, "It seems something went wrong during ASD computation\n");
    fprintf(params->simul_log, "It seems something went wrong during ASD computation\n");
  }
  free(command);
  free(tmp);
}

unsigned **sample_individuals(param *params, arg *args, couple **couples){
  unsigned **wanted = NULL, nb_done = 0;
  int i, j;
  unsigned is_ok;
  char chrom1[SMALL_BUFF_SIZE], chrom2[SMALL_BUFF_SIZE];
  fprintf(stderr, "Sampling individuals for statistics computations...\n");
  fprintf(params->simul_log, "Sampling individuals for statistics computations...\n");
  wanted = allocation_unsigned_matrix(args->sample_size_s1 + args->sample_size_s2 + args->sample_size_adm, 3);
  /* we shuffle chromosome orders for source populations because we used this order to create the couples so we expect a lot of relatedness */
  gsl_ran_shuffle(params->rngs[0], params->order_chroms_s1, params->nb_chrom_s1, sizeof(unsigned));
  gsl_ran_shuffle(params->rngs[0], params->order_chroms_s2, params->nb_chrom_s2, sizeof(unsigned));
  /* first we try to pick unrelated individuals */
  /* we start with admixed population */
  for (i = 0; i < params->current_adm_Ne && nb_done < args->sample_size_adm; i++){
    is_ok = 1;
    for (j = i-1; j < i && j >= 0 && nb_done < args->sample_size_adm; j++){
      if (adm_indivs_share_chromosomes(couples[i], couples[j])){
      	is_ok = 0;
      	break;
      }
      if (is_ok){
	wanted[nb_done][0] = FLAG_ADM;
	wanted[nb_done][1] = 2*i;
	wanted[nb_done][2] = 2*i + 1;
	nb_done ++;
      }
    }
  }
  /* then individuals from s1 population */
  if (nb_done == args->sample_size_adm){
    for (i = 0; i < params->nb_chrom_s1 / 2 && nb_done < args->sample_size_adm + args->sample_size_s1 ; i++){
      /* create chroms identifiers */
      sprintf(chrom1, "s1_%u", params->order_chroms_s1[2*i]);
      sprintf(chrom2, "s1_%u", params->order_chroms_s1[2*i+1]);
      /* check against sampled in admix pop */
      if (adm_shares_chromosome_with_source(chrom1, chrom2, wanted, couples, args->sample_size_adm))
	continue;
      /* if OK, we keep him */
      wanted[nb_done][0] = FLAG_S1;
      wanted[nb_done][1] = params->order_chroms_s1[2*i];
      wanted[nb_done][2] = params->order_chroms_s1[2*i+1];
      nb_done ++;
    }
  }
  /* then individuals from s2 population */
  if (nb_done == args->sample_size_adm + args->sample_size_s1){
    for (i = 0; i < params->nb_chrom_s2 / 2 && nb_done < args->sample_size_adm + args->sample_size_s1 + args->sample_size_s2; i++){
      /* create chroms identifiers */
      sprintf(chrom1, "s2_%u", params->order_chroms_s2[2*i]);
      sprintf(chrom2, "s2_%u", params->order_chroms_s2[2*i+1]);
      /* check against sampled in admix pop */
      if (adm_shares_chromosome_with_source(chrom1, chrom2, wanted, couples, args->sample_size_adm))
	continue;
      /* if OK, we keep him */
      wanted[nb_done][0] = FLAG_S2;
      wanted[nb_done][1] = params->order_chroms_s2[2*i];
      wanted[nb_done][2] = params->order_chroms_s2[2*i+1];
      nb_done ++;
    }
  }
  if (nb_done == args->sample_size_adm + args->sample_size_s1 + args->sample_size_s2)
    return wanted;
  /* if we cannot find unrelated individuals */
  fprintf(stderr, "Cannot find enough unrelated individuals, I will pick samples at random\n");
  fprintf(params->simul_log, "Cannot find enough unrelated individuals, I will pick samples at random\n");
  nb_done = 0;
  for (i = 0; i < args->sample_size_adm; i++){
    wanted[nb_done][0] = FLAG_ADM;
    wanted[nb_done][1] = 2*i;
    wanted[nb_done][2] = 2*i+1;
    nb_done++;
  }
  for (i = 0; i < args->sample_size_s1; i++){
    wanted[nb_done][0] = FLAG_S1;
    wanted[nb_done][1] = params->order_chroms_s1[2*i];
    wanted[nb_done][2] = params->order_chroms_s1[2*i+1];
    nb_done++;
  }
  for (i = 0; i < args->sample_size_s2; i++){
    wanted[nb_done][0] = FLAG_S2;
    wanted[nb_done][1] = params->order_chroms_s2[2*i];
    wanted[nb_done][2] = params->order_chroms_s2[2*i+1];
    nb_done++;
  }
  return wanted;
}

unsigned adm_shares_chromosome_with_source(char *chrom1, char *chrom2, unsigned **wanted, couple **couples, unsigned nb_adm){
  unsigned i, j;
  couple *tmp = NULL;
  for (i = 0; i < nb_adm; i++){
    tmp = couples[wanted[i][1]/2];
    for (j = 0; j < 4; j++){
      if (strncmp(chrom1, tmp->chroms_names[j], SMALL_BUFF_SIZE) == 0 || strncmp(chrom2, tmp->chroms_names[j], SMALL_BUFF_SIZE) == 0){
	return 1;
      }
    }
  }
  return 0;
}

unsigned adm_indivs_share_chromosomes(couple *couple1, couple *couple2){
  unsigned i, j;
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      if (strncmp(couple1->chroms_names[i], couple2->chroms_names[j], SMALL_BUFF_SIZE) == 0)
	return 1;
  return 0;
}

void make_multithreaded_reproduction(param *params, couple **couples, unsigned nb_thread, unsigned nb_snp){
  /* we have to deal with multiple threading here */
  pthread_t *threads = NULL;
  thread_datast *thread_datas = NULL;
  unsigned i, increment, idx_deb, idx_fin;
  int dummy;
  threads = allocation_pthread_t_vector(nb_thread);
  thread_datas = allocation_thread_datast_vector(nb_thread);
  increment = params->current_adm_Ne / nb_thread;
  for (i = 0; i < nb_thread; i++){
    idx_deb = increment * i;
    idx_fin = increment * (i+1) - 1;
    if (i == nb_thread - 1)
      idx_fin = params->current_adm_Ne - 1;
    thread_datas[i].idx_min = idx_deb;
    thread_datas[i].idx_max = idx_fin;
    thread_datas[i].nb_snp = nb_snp;
    thread_datas[i].couples = couples;
    thread_datas[i].params = params;
    thread_datas[i].rng = params->rngs[i];
    dummy = pthread_create(&threads[i], NULL, make_reproduction, (void*) &thread_datas[i]);
    if (dummy != 0)
      exit_on_error(NULL, PTHREAD_ERROR);
  }
  for (i = 0; i < nb_thread; i++){
    pthread_join(threads[i], NULL);
  }
  free(threads);
  free(thread_datas);
}

void *make_reproduction(void *ptr){
  thread_datast *data = (thread_datast*) ptr;
  unsigned idx_couple, idx_snp, real_idx_snp1 = 0, real_idx_snp2 = 0;
  couple *tmp_couple;
  for (idx_couple = data->idx_min; idx_couple <= data->idx_max; idx_couple++){
    tmp_couple = data->couples[idx_couple];
    for (idx_snp = 0; idx_snp < data->nb_snp; idx_snp++){
      real_idx_snp1 = tmp_couple->chroms_names[0][0] == 'a' ? idx_snp : data->params->wanted_snps[idx_snp];
      real_idx_snp2 = tmp_couple->chroms_names[2][0] == 'a' ? idx_snp : data->params->wanted_snps[idx_snp];
      data->params->genos_adm_new[2*idx_couple][idx_snp] = tmp_couple->genos1[(gsl_rng_uniform(data->rng) < .5 ? 0 : 1)][real_idx_snp1];
      data->params->genos_adm_new[2*idx_couple+1][idx_snp] = tmp_couple->genos2[(gsl_rng_uniform(data->rng) < .5 ? 0 : 1)][real_idx_snp2];
    }
  }
  return NULL;
}

void create_one_couple(couple *to_create, couple_param *couple_params, gsl_rng *r){
  unsigned idx1, idx2, idx_chrom11, idx_chrom12, idx_chrom21, idx_chrom22;
  idx1 = gsl_rng_uniform_int(r, couple_params->nb_indiv_pop1);
  idx2 = gsl_rng_uniform_int(r, couple_params->nb_indiv_pop2);
  /* if we draw individuals from the same source population, we do not want an individual to reproduce with himself */
  if (couple_params->genos_pop1 == couple_params->genos_pop2){
    while (idx1 == idx2){
      idx1 = gsl_rng_uniform_int(r, couple_params->nb_indiv_pop1);
      idx2 = gsl_rng_uniform_int(r, couple_params->nb_indiv_pop2);
    }
  }
  idx_chrom11 = couple_params->order_pop1[2*idx1];
  idx_chrom12 = couple_params->order_pop1[2*idx1+1];
  idx_chrom21 = couple_params->order_pop2[2*idx2];
  idx_chrom22 = couple_params->order_pop2[2*idx2+1];
  to_create->genos1[0] = couple_params->genos_pop1[idx_chrom11];
  to_create->genos1[1] = couple_params->genos_pop1[idx_chrom12];
  to_create->genos2[0] = couple_params->genos_pop2[idx_chrom21];
  to_create->genos2[1] = couple_params->genos_pop2[idx_chrom22];
  sprintf(to_create->chroms_names[0], "%s_%i", couple_params->name_pop1, idx_chrom11);
  sprintf(to_create->chroms_names[1], "%s_%i", couple_params->name_pop1, idx_chrom12);
  sprintf(to_create->chroms_names[2], "%s_%i", couple_params->name_pop2, idx_chrom21);
  sprintf(to_create->chroms_names[3], "%s_%i", couple_params->name_pop2, idx_chrom22);
}

void set_couple_param(couple_param *p, unsigned couple_type, param *params){
  switch (couple_type){
  case S1S1 :
    p->genos_pop1 = params->genos_s1;
    p->genos_pop2 = params->genos_s1;
    p->nb_indiv_pop1 = params->nb_chrom_s1 / 2;
    p->nb_indiv_pop2 = params->nb_chrom_s1 / 2;
    p->name_pop1 = strncpy(p->name_pop1, "s1", SMALL_BUFF_SIZE);
    p->name_pop2 = strncpy(p->name_pop2, "s1", SMALL_BUFF_SIZE);
    p->order_pop1 = params->order_chroms_s1;
    p->order_pop2 = params->order_chroms_s1;
    break;
  case S2S2 :
    p->genos_pop1 = params->genos_s2;
    p->genos_pop2 = params->genos_s2;
    p->nb_indiv_pop1 = params->nb_chrom_s2 / 2;
    p->nb_indiv_pop2 = params->nb_chrom_s2 / 2;
    p->name_pop1 = strncpy(p->name_pop1, "s2", SMALL_BUFF_SIZE);
    p->name_pop2 = strncpy(p->name_pop2, "s2", SMALL_BUFF_SIZE);
    p->order_pop1 = params->order_chroms_s2;
    p->order_pop2 = params->order_chroms_s2;
    break;
  case S1S2 :
    p->genos_pop1 = params->genos_s1;
    p->genos_pop2 = params->genos_s2;
    p->nb_indiv_pop1 = params->nb_chrom_s1 / 2;
    p->nb_indiv_pop2 = params->nb_chrom_s2 / 2;
    p->name_pop1 = strncpy(p->name_pop1, "s1", SMALL_BUFF_SIZE);
    p->name_pop2 = strncpy(p->name_pop2, "s2", SMALL_BUFF_SIZE);
    p->order_pop1 = params->order_chroms_s1;
    p->order_pop2 = params->order_chroms_s2;
    break;
  case S1H :
    p->genos_pop1 = params->genos_s1;
    p->genos_pop2 = params->genos_adm_old;
    p->nb_indiv_pop1 = params->nb_chrom_s1 / 2;
    p->nb_indiv_pop2 = params->prev_adm_Ne / 2;
    p->name_pop1 = strncpy(p->name_pop1, "s1", SMALL_BUFF_SIZE);
    p->name_pop2 = strncpy(p->name_pop2, "adm", SMALL_BUFF_SIZE);
    p->order_pop1 = params->order_chroms_s1;
    p->order_pop2 = params->order_chroms_adm;
    break;
  case S2H :
    p->genos_pop1 = params->genos_s2;
    p->genos_pop2 = params->genos_adm_old;
    p->nb_indiv_pop1 = params->nb_chrom_s2 / 2;
    p->nb_indiv_pop2 = params->prev_adm_Ne / 2;
    p->name_pop1 = strncpy(p->name_pop1, "s2", SMALL_BUFF_SIZE);
    p->name_pop2 = strncpy(p->name_pop2, "adm", SMALL_BUFF_SIZE);
    p->order_pop1 = params->order_chroms_s2;
    p->order_pop2 = params->order_chroms_adm;
    break;
  case HH :
    p->genos_pop1 = params->genos_adm_old;
    p->genos_pop2 = params->genos_adm_old;
    p->nb_indiv_pop1 = params->prev_adm_Ne / 2;
    p->nb_indiv_pop2 = params->prev_adm_Ne / 2;
    p->name_pop1 = strncpy(p->name_pop1, "adm", SMALL_BUFF_SIZE);
    p->name_pop2 = strncpy(p->name_pop2, "adm", SMALL_BUFF_SIZE);
    p->order_pop1 = params->order_chroms_adm;
    p->order_pop2 = params->order_chroms_adm;
    break;
  }
}

void reset_couple_param(couple_param *p){
  p->genos_pop1 = NULL;
  p->genos_pop2 = NULL;
  p->nb_indiv_pop1 = 0;
  p->nb_indiv_pop2 = 0;
  p->name_pop1[0] = '\0';
  p->name_pop2[0] = '\0';
  p->order_pop1 = NULL;
  p->order_pop2 = NULL;
}

unsigned *initialize_chrom_order(unsigned nb_chrom){
  unsigned *order = NULL;
  unsigned i;
  order = allocation_unsigned_vector(nb_chrom);
  for (i = 0; i < nb_chrom; i++)
    order[i] = i;
  return order;
}

void count_nb_couple_each_class(unsigned *nb_couples, unsigned Ne, double p1, double p2, gsl_rng *rng){
  double adm, probas[6];
  unsigned i, s = 0, to_change;
  adm = 1.0 - p1 - p2;
  probas[S1S1] = p1 * p1;
  probas[S1H] = 2 * p1 * adm;
  probas[HH] = adm * adm;
  probas[S2H] = 2 * p2 * adm;
  probas[S2S2] = p2 * p2;
  probas[S1S2] = 2 * p1 * p2;
  
  for (i = 0; i < 6; i++){
    nb_couples[i] = (unsigned) round(Ne * probas[i]);
    s += nb_couples[i];
    
  }
  /* correct problems caused by rounding errors */
  /* we add or remove one couple from a random type of couples with non zero probability of occuring */
  while (s != Ne){
    to_change = gsl_rng_uniform_int(rng, 6);
    if (probas[to_change] < 1e-8)
      continue;
    if (s > Ne && nb_couples[to_change] > 0)
      (nb_couples[to_change])--;
    else if (s < Ne)
      (nb_couples[to_change])++;
    s = 0;
    for (i = 0; i < 6; i++)
      s += nb_couples[i];
  }
}
