#include "misc.h"
#include "sumstats.h"
#include "io.h"
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rstat.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_interp.h>
#include <math.h>

/* computes all sumstats... */
void compute_all_sumstats(param *params, arg *args, unsigned **wanted_samples){
  char **genotypes = NULL;
  double **asd = NULL, **mds = NULL, *admixture_props = NULL, *admixture_angles = NULL;
  unsigned i, nb_indiv = args->sample_size_s1 + args->sample_size_s2 + args->sample_size_adm;
  sumstats_datast *sumstats = NULL;
  freqs_datast *frequencies = NULL;
  sumstats = (sumstats_datast*) malloc(sizeof(sumstats_datast));
  initialize_sumstats(sumstats);
  check_alloc(sumstats);
  genotypes = simplify_genotypes(params, args, wanted_samples);
  asd = compute_asd_matrix_local(params, args, genotypes);

  /* stats derived directly from ASD matrix */
  /* they need to be computed NOW because the MDS computation destroys some of the values */
  compute_ASD_stats(asd, args, sumstats);

  
  mds = compute_mds(asd, nb_indiv);
  admixture_props = compute_admixture_proportions(mds, args, nb_indiv);
  admixture_angles = compute_admixture_angles(mds, args, nb_indiv);
  
  /* moments of admixture proportions */
  sumstats->adm_prop_mean = gsl_stats_mean(admixture_props, 1, args->sample_size_adm);
  sumstats->adm_prop_sd = gsl_stats_sd_m(admixture_props, 1, args->sample_size_adm, sumstats->adm_prop_mean);
  sumstats->adm_prop_skew = gsl_stats_skew_m_sd(admixture_props, 1, args->sample_size_adm, sumstats->adm_prop_mean, sumstats->adm_prop_sd);
  sumstats->adm_prop_kurt = gsl_stats_kurtosis_m_sd(admixture_props, 1, args->sample_size_adm, sumstats->adm_prop_mean, sumstats->adm_prop_sd);
  /* for other stats we need to sort admixture_props */
  gsl_sort(admixture_props, 1, args->sample_size_adm);
  sumstats->adm_prop_mode = compute_distribution_mode(admixture_props, args->sample_size_adm);
  
  /* percentiles of admixture_props */
  for (i = 0; i < NB_PERCENTILE; i++)
    sumstats->percentiles[i] = gsl_stats_quantile_from_sorted_data(admixture_props, 1, args->sample_size_adm, i/(NB_PERCENTILE-1.0));

  /* moments of admixture angles */
  sumstats->adm_ang_mean = gsl_stats_mean(admixture_angles, 1, args->sample_size_adm);
  sumstats->adm_ang_sd = gsl_stats_sd_m(admixture_angles, 1, args->sample_size_adm, sumstats->adm_ang_mean);
  sumstats->adm_ang_skew = gsl_stats_skew_m_sd(admixture_angles, 1, args->sample_size_adm, sumstats->adm_ang_mean, sumstats->adm_ang_sd);
  sumstats->adm_ang_kurt = gsl_stats_kurtosis_m_sd(admixture_angles, 1, args->sample_size_adm, sumstats->adm_ang_mean, sumstats->adm_ang_sd);
  /* for other stats we need to sort admixture_props */
  gsl_sort(admixture_angles, 1, args->sample_size_adm);
  sumstats->adm_ang_mode = compute_distribution_mode(admixture_angles, args->sample_size_adm);

  /* percentiles of admixture_angles */
  for (i = 0; i < NB_PERCENTILE; i++)
    sumstats->ang_percentiles[i] = gsl_stats_quantile_from_sorted_data(admixture_angles, 1, args->sample_size_adm, i/(NB_PERCENTILE-1.0));


  /* stats based on allelic frequencies */
  frequencies = compute_allelic_frequencies_local(genotypes, args, sumstats);

  /* heterozygosity statistics */
  sumstats->mean_het_adm = gsl_stats_mean(frequencies->hets_adm, 1, args->nb_snp);
  sumstats->mean_het_s1 = gsl_stats_mean(frequencies->hets_s1, 1, args->nb_snp);
  sumstats->mean_het_s2 = gsl_stats_mean(frequencies->hets_s2, 1, args->nb_snp);
  sumstats->var_het_adm = gsl_stats_variance_m(frequencies->hets_adm, 1, args->nb_snp, sumstats->mean_het_adm);
  sumstats->var_het_s1 = gsl_stats_variance_m(frequencies->hets_s1, 1, args->nb_snp, sumstats->mean_het_s1);
  sumstats->var_het_s2 = gsl_stats_variance_m(frequencies->hets_s2, 1, args->nb_snp, sumstats->mean_het_s2);

  print_sumstats(sumstats, args, params);
  
  free_freqs_datast(frequencies);
  free_char_matrix(genotypes, nb_indiv);
  genotypes = NULL;
  free_double_matrix(asd, nb_indiv);
  asd = NULL;
  free_double_matrix(mds, nb_indiv);
  mds = NULL;
  free(admixture_props);
  free(admixture_angles);
  free(sumstats->percentiles);
  free(sumstats->ang_percentiles);
  free(sumstats);
}

/* computes allelic frequencies, and other things based on genotypes which are useful */
freqs_datast *compute_allelic_frequencies_local(char **genotypes, arg *args, sumstats_datast *sumstats){
  freqs_datast *freqs = NULL;
  thread_datast_freqs *threads_data = NULL;
  pthread_t *threads = NULL;
  unsigned i;
  int dummy;
  double num_f3 = 0, denom_f3 = 0, num_fst_adm_s1 = 0, denom_fst_adm_s1 = 0, num_fst_adm_s2 = 0, denom_fst_adm_s2 = 0, num_fst_s1_s2 = 0, denom_fst_s1_s2 = 0;
  threads_data = (thread_datast_freqs*) malloc(args->nb_thread * sizeof(thread_datast_freqs));
  check_alloc(threads_data);
  threads = allocation_pthread_t_vector(args->nb_thread);
  freqs = allocation_freqs_datast(args);

  for (i = 0; i < args->nb_thread; i++){
    /* could be optimized to deal with nb_snp % nb_thread... */
    threads_data[i].idx_deb = i * args->nb_snp / args->nb_thread;
    threads_data[i].idx_fin = (i+1) * args->nb_snp / args->nb_thread - 1;
    if (i == args->nb_thread - 1)
      threads_data[i].idx_fin = args->nb_snp - 1;
    threads_data[i].freqs = freqs;
    threads_data[i].genotypes = genotypes;
    threads_data[i].args = args;
    initialize_thread_datast_freqs(&threads_data[i], args);

    dummy = pthread_create(&threads[i], NULL, compute_allelic_frequencies_snp_subset, (void*) &threads_data[i]);
    if (dummy != 0)
      exit_on_error(NULL, PTHREAD_ERROR);
  }
  for (i = 0; i < args->nb_thread; i++)
    pthread_join(threads[i], NULL);

  /* finalizing F3 */
  for (i = 0; i < args->nb_thread; i++){
    num_f3 += threads_data[i].accu_num_f3;
    denom_f3 += threads_data[i].accu_denom_f3;
  }
  sumstats->f3 = num_f3 / denom_f3;

  /* finalizing Fst */
  for (i = 0; i < args->nb_thread; i++){
    num_fst_adm_s1 += threads_data[i].accu_num_fst_adm_s1;
    denom_fst_adm_s1 += threads_data[i].accu_denom_fst_adm_s1;
    num_fst_adm_s2 += threads_data[i].accu_num_fst_adm_s2;
    denom_fst_adm_s2 += threads_data[i].accu_denom_fst_adm_s2;
    num_fst_s1_s2 += threads_data[i].accu_num_fst_s1_s2;
    denom_fst_s1_s2 += threads_data[i].accu_denom_fst_s1_s2;
  }
  sumstats->fst_adm_s1 = num_fst_adm_s1 / denom_fst_adm_s1;
  sumstats->fst_adm_s2 = num_fst_adm_s2 / denom_fst_adm_s2;
  sumstats->fst_s1_s2 = num_fst_s1_s2 / denom_fst_s1_s2;

  finalize_inbreeding_computations(threads_data, sumstats, args);
  
  free(threads);
  for (i = 0; i < args->nb_thread; i++){
    free(threads_data[i].nb_het_sites_adm);
    free(threads_data[i].nb_het_sites_s1);
    free(threads_data[i].nb_het_sites_s2);
  }
  free(threads_data);
  return freqs;
}

void initialize_thread_datast_freqs(thread_datast_freqs *data, arg *args){
  unsigned i;
  double tmp1, tmp2;
  data->accu_num_f3 = 0;
  data->accu_denom_f3 = 0;
  data->nb_sites_considered_F_adm = 0;
  data->nb_sites_considered_F_s1 = 0;
  data->nb_sites_considered_F_s2 = 0;
  data->nb_het_sites_adm = allocation_unsigned_vector(args->sample_size_adm);
  data->nb_het_sites_s1 = allocation_unsigned_vector(args->sample_size_s1);
  data->nb_het_sites_s2 = allocation_unsigned_vector(args->sample_size_s2);
  for (i = 0; i < args->sample_size_adm; i++)
    data->nb_het_sites_adm[i] = 0;
  for (i = 0; i < args->sample_size_s1; i++)
    data->nb_het_sites_s1[i] = 0;
  for (i = 0; i < args->sample_size_s2; i++)
    data->nb_het_sites_s2[i] = 0;
  data->expected_hom_adm = 0;
  data->expected_hom_s1 = 0;
  data->expected_hom_s2 = 0;
  data->n_bar_adm_s1 = 0.5 * (args->sample_size_adm + args->sample_size_s1);
  data->n_bar_adm_s2 = 0.5 * (args->sample_size_adm + args->sample_size_s2);
  data->n_bar_s1_s2 = 0.5 * (args->sample_size_s1 + args->sample_size_s2);
  tmp1 = 2 * data->n_bar_adm_s1;
  tmp2 = (args->sample_size_adm * args->sample_size_adm) + (args->sample_size_s1 * args->sample_size_s1);
  data->n_C_adm_s1 = tmp1 - tmp2/tmp1;
  tmp1 = 2 * data->n_bar_adm_s2;
  tmp2 = (args->sample_size_adm * args->sample_size_adm) + (args->sample_size_s2 * args->sample_size_s2);
  data->n_C_adm_s2 = tmp1 - tmp2/tmp1;
  tmp1 = 2 * data->n_bar_s1_s2;
  tmp2 = (args->sample_size_s1 * args->sample_size_s1) + (args->sample_size_s2 * args->sample_size_s2);
  data->n_C_s1_s2 = tmp1 - tmp2/tmp1;
  data->accu_num_fst_adm_s1 = 0;
  data->accu_denom_fst_adm_s1 = 0;
  data->accu_num_fst_adm_s2 = 0;
  data->accu_denom_fst_adm_s2 = 0;
  data->accu_num_fst_s1_s2 = 0;
  data->accu_denom_fst_s1_s2 = 0;

}


/* gives same results as vcftools + R (diff ~ 1e-7) but nothing can be done... */
/* difference comes from the fact that vcftools output precision is 5 digits, which makes a difference when computing mean/variance on small samples */
/* tested by recompiling vcftools with higher output precision -> gives the same result */
void finalize_inbreeding_computations(thread_datast_freqs *threads_data, sumstats_datast *sumstats, arg *args){
  double *F_container = NULL, expected_hom;
  unsigned i, j, nb_site_considered, observed_hom, observed_het;
  /* admixed pop */
  F_container = allocation_double_vector(args->sample_size_adm);
  nb_site_considered = 0;
  expected_hom = 0;
  /* cummulate results from different threads */
  for (j = 0; j < args->nb_thread; j++){
    nb_site_considered += threads_data[j].nb_sites_considered_F_adm;
    expected_hom += threads_data[j].expected_hom_adm;
  }
  for (i = 0; i < args->sample_size_adm; i++){
    observed_het = 0;
    for (j = 0; j < args->nb_thread; j++)
      observed_het += threads_data[j].nb_het_sites_adm[i];
    observed_hom = nb_site_considered - observed_het;
    /* compute F */
    F_container[i] = (observed_hom - expected_hom) / (nb_site_considered - expected_hom);
  }
  /* get the wanted sumstats */
  sumstats->mean_F_adm = gsl_stats_mean(F_container, 1, args->sample_size_adm);
  sumstats->var_F_adm = gsl_stats_variance_m(F_container, 1, args->sample_size_adm, sumstats->mean_F_adm);
  free(F_container);

  /* admixed s1 */
  F_container = allocation_double_vector(args->sample_size_s1);
  nb_site_considered = 0;
  expected_hom = 0;
  for (j = 0; j < args->nb_thread; j++){
    nb_site_considered += threads_data[j].nb_sites_considered_F_s1;
    expected_hom += threads_data[j].expected_hom_s1;
  }
  for (i = 0; i < args->sample_size_s1; i++){
    observed_het = 0;
    for (j = 0; j < args->nb_thread; j++)
      observed_het += threads_data[j].nb_het_sites_s1[i];
    observed_hom = nb_site_considered - observed_het;
    F_container[i] = (observed_hom - expected_hom) / (nb_site_considered - expected_hom);
  }
  sumstats->mean_F_s1 = gsl_stats_mean(F_container, 1, args->sample_size_s1);
  sumstats->var_F_s1 = gsl_stats_variance_m(F_container, 1, args->sample_size_s1, sumstats->mean_F_s1);
  free(F_container);

  /* admixed s2 */
  F_container = allocation_double_vector(args->sample_size_s2);
  nb_site_considered = 0;
  expected_hom = 0;
  for (j = 0; j < args->nb_thread; j++){
    nb_site_considered += threads_data[j].nb_sites_considered_F_s2;
    expected_hom += threads_data[j].expected_hom_s2;
  }
  for (i = 0; i < args->sample_size_s2; i++){
    observed_het = 0;
    for (j = 0; j < args->nb_thread; j++)
      observed_het += threads_data[j].nb_het_sites_s2[i];
    observed_hom = nb_site_considered - observed_het;
    F_container[i] = (observed_hom - expected_hom) / (nb_site_considered - expected_hom);
  }
  sumstats->mean_F_s2 = gsl_stats_mean(F_container, 1, args->sample_size_s2);
  sumstats->var_F_s2 = gsl_stats_variance_m(F_container, 1, args->sample_size_s2, sumstats->mean_F_s2);
  free(F_container);

}

/* computes a lot of things necessary for computation of sumstats based on alleles frequencies */
void *compute_allelic_frequencies_snp_subset(void *ptr){
  thread_datast_freqs *data = (thread_datast_freqs*) ptr;
  double tmp_freq, h_bar_adm_s1, h_bar_adm_s2, h_bar_s1_s2, h_tilde_adm, h_tilde_s1, h_tilde_s2, p_bar_adm_s1, p_bar_adm_s2, p_bar_s1_s2, ssquare_adm_s1, ssquare_adm_s2, ssquare_s1_s2, tmp_a, tmp_b;
  unsigned i, j;

  /* we compute allelic frequency of allele 1 (simplifies code) */
  for (i = data->idx_deb; i <= data->idx_fin; i++){
    tmp_freq = 0;
    h_tilde_adm = 0;
    h_tilde_s1 = 0;
    h_tilde_s2 = 0;
    /* admixed population */
    for (j = 0; j < data->args->sample_size_adm; j++){
      tmp_freq += data->genotypes[j][i]; /* update frequencies */
      if (data->genotypes[j][i] == 1){
	data->nb_het_sites_adm[j]++; /* update number of heterozygous site per individual */
	h_tilde_adm++;		     /* update for Fst */
      }
    }
    h_tilde_adm /= data->args->sample_size_adm;
    data->freqs->freqs_adm[i] = tmp_freq / (2*data->args->sample_size_adm); /* finalize frequency computation */
    data->freqs->hets_adm[i] = 2 * data->freqs->freqs_adm[i] * (1 - data->freqs->freqs_adm[i]) * (2*data->args->sample_size_adm) / ((2*data->args->sample_size_adm) - 1); /* compute expected heterozygosity per site */
    /* update stats to compute inbreeding coefficient */
    if (data->freqs->freqs_adm[i] > 0 && data->freqs->freqs_adm[i] < 1){
      data->nb_sites_considered_F_adm++;
      data->expected_hom_adm += 1 - data->freqs->hets_adm[i];
    }
    /* pop s1 */
    tmp_freq = 0;
    for (j = data->args->sample_size_adm; j < data->args->sample_size_adm + data->args->sample_size_s1; j++){
      tmp_freq += data->genotypes[j][i];
      if (data->genotypes[j][i] == 1){
	data->nb_het_sites_s1[j - data->args->sample_size_adm]++;
	h_tilde_s1++;
      }
    }
    h_tilde_s1 /= data->args->sample_size_s1;
    data->freqs->freqs_s1[i] = tmp_freq / (2*data->args->sample_size_s1);
    data->freqs->hets_s1[i] = 2 * data->freqs->freqs_s1[i] * (1 - data->freqs->freqs_s1[i]) * (2*data->args->sample_size_s1) / (2*data->args->sample_size_s1 - 1);
    if (data->freqs->freqs_s1[i] > 0 && data->freqs->freqs_s1[i] < 1){
      data->nb_sites_considered_F_s1++;
      data->expected_hom_s1 += 1 - data->freqs->hets_s1[i];
    }
    /* pop s2 */
    tmp_freq = 0;
    for (j = data->args->sample_size_adm + data->args->sample_size_s1; j < data->args->sample_size_adm + data->args->sample_size_s1 + data->args->sample_size_s2; j++){
      tmp_freq += data->genotypes[j][i];
      if (data->genotypes[j][i] == 1){
	data->nb_het_sites_s2[j - data->args->sample_size_adm - data->args->sample_size_s1]++;
	h_tilde_s2++;
      }
    }
    h_tilde_s2 /= data->args->sample_size_s2;
    data->freqs->freqs_s2[i] = tmp_freq / (2*data->args->sample_size_s2);
    data->freqs->hets_s2[i] = 2 * data->freqs->freqs_s2[i] * (1 - data->freqs->freqs_s2[i]) * (2*data->args->sample_size_s2) / (2*data->args->sample_size_s2 - 1);
    if (data->freqs->freqs_s2[i] > 0 && data->freqs->freqs_s2[i] < 1){
      data->nb_sites_considered_F_s2++;
      data->expected_hom_s2 += 1 - data->freqs->hets_s2[i];
    }

    /* update accus for f3 */
    data->accu_num_f3 += (data->freqs->freqs_adm[i] - data->freqs->freqs_s1[i]) * (data->freqs->freqs_adm[i] - data->freqs->freqs_s2[i]);
    data->accu_denom_f3 += 2 * data->freqs->freqs_adm[i] * (1 - data->freqs->freqs_adm[i]);

    /* update accus for Fst */
    h_bar_adm_s1 = compute_bar_stat_Fst(data->args->sample_size_adm, data->args->sample_size_s1, h_tilde_adm, h_tilde_s1, data->n_bar_adm_s1);
    h_bar_adm_s2 = compute_bar_stat_Fst(data->args->sample_size_adm, data->args->sample_size_s2, h_tilde_adm, h_tilde_s2, data->n_bar_adm_s2);
    h_bar_s1_s2 = compute_bar_stat_Fst(data->args->sample_size_s1, data->args->sample_size_s2, h_tilde_s1, h_tilde_s2, data->n_bar_s1_s2);
    p_bar_adm_s1 = compute_bar_stat_Fst(data->args->sample_size_adm, data->args->sample_size_s1, data->freqs->freqs_adm[i], data->freqs->freqs_s1[i], data->n_bar_adm_s1);
    p_bar_adm_s2 = compute_bar_stat_Fst(data->args->sample_size_adm, data->args->sample_size_s2, data->freqs->freqs_adm[i], data->freqs->freqs_s2[i], data->n_bar_adm_s2);
    p_bar_s1_s2 = compute_bar_stat_Fst(data->args->sample_size_s1, data->args->sample_size_s2, data->freqs->freqs_s1[i], data->freqs->freqs_s2[i], data->n_bar_s1_s2);
    ssquare_adm_s1 = compute_ssquare_Fst(data->args->sample_size_adm, data->args->sample_size_s1, data->freqs->freqs_adm[i], data->freqs->freqs_s1[i], data->n_bar_adm_s1, p_bar_adm_s1);
    ssquare_adm_s2 = compute_ssquare_Fst(data->args->sample_size_adm, data->args->sample_size_s2, data->freqs->freqs_adm[i], data->freqs->freqs_s2[i], data->n_bar_adm_s2, p_bar_adm_s2);
    ssquare_s1_s2 = compute_ssquare_Fst(data->args->sample_size_s1, data->args->sample_size_s2, data->freqs->freqs_s1[i], data->freqs->freqs_s2[i], data->n_bar_s1_s2, p_bar_s1_s2);

    tmp_a = compute_a_Fst(data->n_bar_adm_s1, data->n_C_adm_s1, ssquare_adm_s1, p_bar_adm_s1, h_bar_adm_s1);
    tmp_b = compute_b_Fst(data->n_bar_adm_s1, ssquare_adm_s1, p_bar_adm_s1, h_bar_adm_s1);
    data->accu_num_fst_adm_s1 += tmp_a;
    data->accu_denom_fst_adm_s1 += tmp_a + tmp_b + h_bar_adm_s1 / 2;

    tmp_a = compute_a_Fst(data->n_bar_adm_s2, data->n_C_adm_s2, ssquare_adm_s2, p_bar_adm_s2, h_bar_adm_s2);
    tmp_b = compute_b_Fst(data->n_bar_adm_s2, ssquare_adm_s2, p_bar_adm_s2, h_bar_adm_s2);
    data->accu_num_fst_adm_s2 += tmp_a;
    data->accu_denom_fst_adm_s2 += tmp_a + tmp_b + h_bar_adm_s2 / 2;

    tmp_a = compute_a_Fst(data->n_bar_s1_s2, data->n_C_s1_s2, ssquare_s1_s2, p_bar_s1_s2, h_bar_s1_s2);
    tmp_b = compute_b_Fst(data->n_bar_s1_s2, ssquare_s1_s2, p_bar_s1_s2, h_bar_s1_s2);
    data->accu_num_fst_s1_s2 += tmp_a;
    data->accu_denom_fst_s1_s2 += tmp_a + tmp_b + h_bar_s1_s2 / 2;

    
  }
  return NULL;
}

/* computes a term in Fst */
double compute_a_Fst(double n_bar, double n_C, double ssquare, double p_bar, double h_bar){
  double a;
  a = p_bar * (1 - p_bar) - ssquare / 2 - h_bar /4;
  a = ssquare - a / (n_bar - 1);
  a = n_bar * a / n_C;
  return a;
}

/* computes b term in Fst */
double compute_b_Fst(double n_bar, double ssquare, double p_bar, double h_bar){
  double b;
  b = p_bar * (1 - p_bar) - ssquare / 2;
  b -= (2 * n_bar - 1) * h_bar / (4 * n_bar);
  b = n_bar * b / (n_bar - 1);
  return b;
}

/* compute variance in allelic frequencies */
double compute_ssquare_Fst(unsigned n1, unsigned n2, double p1, double p2, double n_bar, double p_bar){
  double ssquare;
  ssquare = (n1 * (p1 - p_bar) * (p1 - p_bar) + n2 * (p2 - p_bar) * (p2 - p_bar)) / n_bar;
  return ssquare;
}


/* computes Weir bar stats */
double compute_bar_stat_Fst(unsigned ssize1, unsigned ssize2, double tilde1, double tilde2, double n_bar){
  double stat;
  stat = (ssize1 * tilde1 + ssize2 * tilde2) / (2.0 * n_bar);
  return stat;
}

void free_freqs_datast(freqs_datast *f){
  free(f->freqs_adm);
  free(f->freqs_s1);
  free(f->freqs_s2);
  free(f->hets_adm);
  free(f->hets_s1);
  free(f->hets_s2);
  free(f);
}

freqs_datast *allocation_freqs_datast(arg *args){
  freqs_datast *f = NULL;
  f = (freqs_datast*) malloc(sizeof(freqs_datast));
  check_alloc(f);
  f->freqs_adm = allocation_double_vector(args->nb_snp);
  f->freqs_s1 = allocation_double_vector(args->nb_snp);
  f->freqs_s2 = allocation_double_vector(args->nb_snp);
  f->hets_adm = allocation_double_vector(args->nb_snp);
  f->hets_s1 = allocation_double_vector(args->nb_snp);
  f->hets_s2 = allocation_double_vector(args->nb_snp);
  return f;
}

void initialize_sumstats(sumstats_datast *s){
  unsigned i;
  s->adm_prop_mean = GSL_NAN;
  s->adm_prop_sd = GSL_NAN;
  s->adm_prop_skew = GSL_NAN;
  s->adm_prop_kurt = GSL_NAN;
  s->adm_prop_mode = GSL_NAN;
  s->adm_ang_mean = GSL_NAN;
  s->adm_ang_sd = GSL_NAN;
  s->adm_ang_skew = GSL_NAN;
  s->adm_ang_kurt = GSL_NAN;
  s->adm_ang_mode = GSL_NAN;
  s->mean_asd_adm = GSL_NAN;
  s->var_asd_adm = GSL_NAN;
  s->mean_asd_s1 = GSL_NAN;
  s->var_asd_s1 = GSL_NAN;
  s->mean_asd_s2 = GSL_NAN;
  s->var_asd_s2 = GSL_NAN;
  s->mean_asd_adm_s1 = GSL_NAN;
  s->var_asd_adm_s1 = GSL_NAN;
  s->mean_asd_adm_s2 = GSL_NAN;
  s->var_asd_adm_s2 = GSL_NAN;
  s->mean_asd_s1_s2 = GSL_NAN;
  s->var_asd_s1_s2 = GSL_NAN;
  s->percentiles = allocation_double_vector(NB_PERCENTILE);
  s->ang_percentiles = allocation_double_vector(NB_PERCENTILE);
  for (i = 0; i < NB_PERCENTILE; i++){
    s->percentiles[i] = GSL_NAN;
    s->ang_percentiles[i] = GSL_NAN;
  }
  s->f3 = GSL_NAN;
  s->mean_het_adm = GSL_NAN;
  s->var_het_adm = GSL_NAN;
  s->mean_het_s1 = GSL_NAN;
  s->var_het_s1 = GSL_NAN;
  s->mean_het_s2 = GSL_NAN;
  s->var_het_s2 = GSL_NAN;
  s->mean_F_adm = GSL_NAN;
  s->var_F_adm = GSL_NAN;
  s->mean_F_s1 = GSL_NAN;
  s->var_F_s1 = GSL_NAN;
  s->mean_F_s2 = GSL_NAN;
  s->var_F_s2 = GSL_NAN;
  s->fst_adm_s1 = GSL_NAN;
  s->fst_adm_s2 = GSL_NAN;
  s->fst_s1_s2 = GSL_NAN;
}

/* uses running statistics to prevent too much alloc/free of possibly different sizes */
void compute_ASD_stats(double **asd, arg *args, sumstats_datast *sumstats){
  unsigned i, j, beg_adm, end_adm, beg_s1, end_s1, beg_s2, end_s2;
  gsl_rstat_workspace *w = gsl_rstat_alloc();
  gsl_rstat_reset(w);
  beg_adm = 0;
  end_adm = args->sample_size_adm;
  beg_s1 = end_adm;
  end_s1 = beg_s1 + args->sample_size_s1;
  beg_s2 = end_s1;
  end_s2 = beg_s2 + args->sample_size_s2;
  /* within pop statistics */
  for (i = beg_adm; i < end_adm; i++)
    for (j = i+1; j < end_adm; j++)
      gsl_rstat_add(asd[i][j], w);
  sumstats->mean_asd_adm = gsl_rstat_mean(w);
  sumstats->var_asd_adm = gsl_rstat_variance(w);
  gsl_rstat_reset(w);
  for (i = beg_s1; i < end_s1; i++)
    for (j = i+1; j < end_s1; j++)
      gsl_rstat_add(asd[i][j], w);
  sumstats->mean_asd_s1 = gsl_rstat_mean(w);
  sumstats->var_asd_s1 = gsl_rstat_variance(w);
  gsl_rstat_reset(w);

  for (i = beg_s2; i < end_s2; i++)
    for (j = i+1; j < end_s2; j++)
      gsl_rstat_add(asd[i][j], w);
  sumstats->mean_asd_s2 = gsl_rstat_mean(w);
  sumstats->var_asd_s2 = gsl_rstat_variance(w);
  gsl_rstat_reset(w);

  /* between pop statistics */
  for (i = beg_adm; i < end_adm; i++)
    for (j = beg_s1; j < end_s1; j++)
      gsl_rstat_add(asd[i][j], w);
  sumstats->mean_asd_adm_s1 = gsl_rstat_mean(w);
  sumstats->var_asd_adm_s1 = gsl_rstat_variance(w);
  gsl_rstat_reset(w);
  
  for (i = beg_adm; i < end_adm; i++)
    for (j = beg_s2; j < end_s2; j++)
      gsl_rstat_add(asd[i][j], w);
  sumstats->mean_asd_adm_s2 = gsl_rstat_mean(w);
  sumstats->var_asd_adm_s2 = gsl_rstat_variance(w);
  gsl_rstat_reset(w);

  for (i = beg_s1; i < end_s1; i++)
    for (j = beg_s2; j < end_s2; j++)
      gsl_rstat_add(asd[i][j], w);
  sumstats->mean_asd_s1_s2 = gsl_rstat_mean(w);
  sumstats->var_asd_s1_s2 = gsl_rstat_variance(w);
  gsl_rstat_reset(w);

  gsl_rstat_free(w);
}

/* computes the KDE of a distribution and returns the mode */
/* gives the same result as R (diff ~ 1e-9) */
double compute_distribution_mode(double *data, unsigned size){
  double best_x = -1000, best_y = -1, **kde;
  unsigned i;

  kde = compute_kernel_density_estimate(data, size);

  for (i = 0; i < NB_POINTS_KDE; i++){
    if (kde[i][1] > best_y){
      best_y = kde[i][1];
      best_x = kde[i][0];
    }
  }
  free_double_matrix(kde, NB_POINTS_KDE);
  return best_x;
}

/* computes the KDE the same way as R, using R default */
/* could be much simpler (i.e. sum of gaussians), but needs to give the same result as R because it's what we are using on the real data... */
/* this gives the same result as R, up to machine epsilon (1e-9) */
double **compute_kernel_density_estimate(double *data, unsigned size){
  double bw, from, to, lo, up, *y = NULL, *kords = NULL, new_real, new_imag, *xords = NULL, delta, tmp_x, **result = NULL;
  unsigned i;
  gsl_interp *workspace = NULL;
  gsl_interp_accel *acc = NULL;
  
  /* compute bandwidth and ranges limits */
  bw = nrd0(data, size);
  from = data[0] - 3*bw;
  to = data[size-1] + 3*bw;
  lo = from - 4*bw;
  up = to + 4*bw;
  
  y = compute_bin_dist(data, lo, up, size);
  kords = compute_kords(up, lo, bw);

  gsl_fft_real_radix2_transform(y, 1, 2*NB_POINTS_KDE);
  gsl_fft_real_radix2_transform(kords, 1, 2*NB_POINTS_KDE);
  
  /* conjugate of fft(kords) */
  for (i = NB_POINTS_KDE+1; i < 2*NB_POINTS_KDE; i++)
    kords[i] = -kords[i];
  
  /* fft(y) * Conj(fft(kords)) -> store in kords */
  for (i = 0; i <= NB_POINTS_KDE; i++){
    /* real part */
    new_real = kords[i]*y[i];
    if (i != 0 && i != NB_POINTS_KDE){
      new_real -= kords[2*NB_POINTS_KDE-i] * y[2*NB_POINTS_KDE-i];
      /* imaginary part */
      new_imag = kords[i] * y[2*NB_POINTS_KDE-i] + y[i] * kords[2*NB_POINTS_KDE-i];
    }
    kords[i] = new_real;
    if (i != 0 && i != NB_POINTS_KDE){
      kords[2*NB_POINTS_KDE-i] = new_imag;
    }
  }
  
  /* inverse fft */
  gsl_fft_halfcomplex_radix2_inverse(kords, 1, 2*NB_POINTS_KDE);

  /* wanted data are in 0..NB_POINTS_KDE and should be >= 0 */
  for (i = 0; i < NB_POINTS_KDE; i++)
    kords[i] = GSL_MAX(kords[i], 0);

  /* finally, interpolation */
  xords = allocation_double_vector(NB_POINTS_KDE);
  tmp_x = lo;
  delta = (up - lo) / (NB_POINTS_KDE-1);
  for (i = 0; i < NB_POINTS_KDE; i++){
    xords[i] = tmp_x;
    tmp_x += delta;
  }
  
  result = allocation_double_matrix(NB_POINTS_KDE, 2);
  tmp_x = from;
  delta = (to-from) / (NB_POINTS_KDE-1);
  for (i = 0; i < NB_POINTS_KDE; i++){
    result[i][0] = tmp_x;
    tmp_x += delta;
  }
  /* we use linear interpolation (same as R) */
  workspace = gsl_interp_alloc(gsl_interp_linear, NB_POINTS_KDE);
  gsl_interp_init(workspace, xords, kords, NB_POINTS_KDE);
  acc = gsl_interp_accel_alloc();
  gsl_interp_accel_reset(acc);

  for (i = 0; i < NB_POINTS_KDE; i++)
    result[i][1] = gsl_interp_eval(workspace, xords, kords, result[i][0], acc);
  
  gsl_interp_free(workspace);
  gsl_interp_accel_free(acc);
  free(y);
  free(kords);
  free(xords);
  return result;
}

double *compute_kords(double hi, double lo, double bw){
  double *kords = NULL, delta, tmp = 0;
  unsigned i;
  kords = allocation_double_vector(2*NB_POINTS_KDE);
  delta = 2*(hi-lo) / (2*NB_POINTS_KDE-1);
  for (i = 0; i < 2*NB_POINTS_KDE; i++){
    kords[i] = tmp;
    tmp += delta;
  }
  for (i = 1; i < NB_POINTS_KDE; i++)
    kords[2*NB_POINTS_KDE-i] = -kords[i];
  for (i = 0; i < 2*NB_POINTS_KDE; i++)
    kords[i] = gsl_ran_gaussian_pdf(kords[i], bw);
  return kords;
}

/* gives same result as R, up to machine eps (~1e-9) */
double *compute_bin_dist(double *data, double lo, double hi, unsigned size){
  double *y = NULL, xdelta, xpos, fx, weight;
  unsigned i, ixmin, ixmax;
  int ix;
  y = allocation_double_vector(2*NB_POINTS_KDE);
  for (i = 0; i < 2*NB_POINTS_KDE; i++)
    y[i] = 0;
  ixmin = 0;
  ixmax = NB_POINTS_KDE - 2;
  xdelta = (hi - lo) / (NB_POINTS_KDE - 1);
  weight = 1.0/size;
  
  for (i = 0; i < size; i++){
    xpos = (data[i] - lo) / xdelta;
    ix = (int) floor(xpos);
    fx = xpos - ix;
    if (ixmin <= ix && ix <= ixmax){
      y[ix] += (1 - fx) * weight;
      y[ix+1] += fx * weight;
    }
    else if (ix == -1)
      y[0] += fx * weight;
    else if(ix == ixmax+1)
      y[ix] += (1 - fx) * weight;
  }
  
  return y;
}

/* computes bandwidth for KDE */
double nrd0(double *data, unsigned size){
  double hi, iqr, lo;
  hi = gsl_stats_sd(data, 1, size);
  iqr = gsl_stats_quantile_from_sorted_data(data, 1, size, 0.75) - gsl_stats_quantile_from_sorted_data(data, 1, size, 0.25);
  lo = GSL_MIN(hi, iqr/1.34);
  return 0.9 * lo * pow(size, -0.2);
}

/* computes angles for each adm indiv on MDS projection */
double *compute_admixture_angles(double **mds, arg *args, unsigned nb_indiv){
  double *angles = NULL, *centroid_s1 = NULL, *centroid_s2 = NULL, *vec_a = NULL, *vec_b = NULL;
  double norm_vec_a, norm_vec_b, cosine;
  unsigned i, j;
  fprintf(stderr, "Computing admixture angles...\n");
  vec_a = allocation_double_vector(NB_DIM_MDS);
  vec_b = allocation_double_vector(NB_DIM_MDS);
  /* the order of individuals in the MDS is the same as in the genotypes:
     adm...s1...s2 */
  /* compute centroids */
  centroid_s1 = compute_centroid(mds, args->sample_size_adm, args->sample_size_adm + args->sample_size_s1 - 1);
  centroid_s2 = compute_centroid(mds, args->sample_size_adm + args->sample_size_s1, nb_indiv-1);

  angles = allocation_double_vector(args->sample_size_adm);
  /* loop over adm indiv */
  for (i = 0; i < args->sample_size_adm; i++){
    /* compute vector to both centroids */
    for (j = 0; j < NB_DIM_MDS; j++){
      vec_a[j] = mds[i][j] - centroid_s1[j];
      vec_b[j] = mds[i][j] - centroid_s2[j];
    }
    norm_vec_a = sqrt(compute_dot_product(vec_a, vec_a));
    norm_vec_b = sqrt(compute_dot_product(vec_b, vec_b));
    cosine = compute_dot_product(vec_a, vec_b) / (norm_vec_a * norm_vec_b);
    angles[i] = acos(cosine);
  }
  
  free(centroid_s1);
  free(centroid_s2);
  free(vec_a);
  free(vec_b);
  return angles;
}

/* computes admixture proportions based on MDS */
double *compute_admixture_proportions(double **mds, arg *args, unsigned nb_indiv){
  double *props = NULL, *centroid_s1 = NULL, *centroid_s2 = NULL, centroid_dist, *vec_b = NULL, *vec_a = NULL, *proj = NULL, tmp, dot_prod_vec_b_vec_b;
  unsigned i, j;
  fprintf(stderr, "Computing admixture proportions...\n");
  /* the order of individuals in the MDS is the same as in the genotypes:
     adm...s1...s2 */
  /* compute centroids */
  centroid_s1 = compute_centroid(mds, args->sample_size_adm, args->sample_size_adm + args->sample_size_s1 - 1);
  centroid_s2 = compute_centroid(mds, args->sample_size_adm + args->sample_size_s1, nb_indiv-1);

  /* distance between centroids */
  centroid_dist = compute_dist_two_points(centroid_s1, centroid_s2);

  /* vector between centroids */
  vec_b = compute_vector_two_points(centroid_s1, centroid_s2);
  dot_prod_vec_b_vec_b = compute_dot_product(vec_b, vec_b);
  
  props = allocation_double_vector(args->sample_size_adm);
  /* loop over admixted individuals */
  vec_a = allocation_double_vector(NB_DIM_MDS);
  proj = allocation_double_vector(NB_DIM_MDS);
  for (i = 0; i < args->sample_size_adm; i++){
    /* vector from centroid_s1 to indiv[i] */
    for (j = 0; j < NB_DIM_MDS; j++)
      vec_a[j] = mds[i][j] - centroid_s1[j];
    /* project vec_a on vec_b */
    tmp = compute_dot_product(vec_a, vec_b) / dot_prod_vec_b_vec_b;
    for (j = 0; j < NB_DIM_MDS; j++)
      proj[j] = centroid_s1[j] + tmp * vec_b[j];
    /* compute proportion */
    props[i] = 1 - compute_dist_two_points(proj, centroid_s1) / centroid_dist;
  }
  
  free(vec_a);  
  free(centroid_s1);
  free(centroid_s2);
  free(vec_b);
  free(proj);
  return props;
}


double compute_dot_product(double *v1, double *v2){
  double dot = 0.0;
  unsigned i;
  for (i = 0; i < NB_DIM_MDS; i++)
    dot += v1[i] * v2[i];
  return dot;
}

/* computes the vector between two points */
double *compute_vector_two_points(double *source, double *dest){
  double *v = NULL;
  unsigned i;
  v = allocation_double_vector(NB_DIM_MDS);
  for (i = 0; i < NB_DIM_MDS; i++)
    v[i] = dest[i] - source[i];
  return v;
}

/* computes the distance between two points */
double compute_dist_two_points(double *p1, double *p2){
  double dist = 0, tmp;
  unsigned i;
  for (i = 0; i < NB_DIM_MDS; i++){
    tmp = p2[i] - p1[i];
    dist += tmp * tmp;
  }
  return sqrt(dist);
}

/* computes coordinates of centroid for a set of points */
double *compute_centroid(double **coords, unsigned idx_min, unsigned idx_max){
  double *centroid = NULL, tmp;
  unsigned i, j;
  centroid = allocation_double_vector(NB_DIM_MDS);
  for (j = 0; j < NB_DIM_MDS; j++){
    tmp = 0.0;
    for (i = idx_min; i <= idx_max; i++)
      tmp += coords[i][j];
    centroid[j] = tmp / (idx_max - idx_min + 1);
  }
  
  return centroid;
}

/* computes a classical MDS */
/* gives the same results as R "cmdscale" function (up to the sign, which is normal) */
double **compute_mds(double **dists, unsigned nb_indiv){
  gsl_matrix *M = NULL;
  size_t i, j;
  double tmp;
  gsl_vector *eval = NULL;
  gsl_matrix *evec = NULL;
  gsl_eigen_symmv_workspace *w = NULL;
  double **mds;
  
  fprintf(stderr, "Computing MDS...\n");
  
  M = gsl_matrix_alloc(nb_indiv, nb_indiv);

  /* set values in M, M=dists^2 */
  for (i = 0; i < nb_indiv; i++){
    gsl_matrix_set(M, i, i, 0.0);
    for (j = i+1; j < nb_indiv; j++){
      tmp = dists[i][j] * dists[i][j];
      gsl_matrix_set(M, i, j, tmp);
      gsl_matrix_set(M, j, i, tmp);
    }
  }

  /* double centering */
  for (i = 0; i < nb_indiv; i++){
    tmp = 0;
    for (j = 0; j < nb_indiv; j++)
      tmp += gsl_matrix_get(M, i, j);
    tmp /= nb_indiv;
    for (j = 0; j < nb_indiv; j++)
      gsl_matrix_set(M, i, j, gsl_matrix_get(M, i, j) - tmp);
  }
  for (j = 0; j < nb_indiv; j++){
    tmp = 0;
    for (i = 0; i < nb_indiv; i++)
      tmp += gsl_matrix_get(M, i, j);
    tmp /= nb_indiv;
    for (i = 0; i < nb_indiv; i++)
      gsl_matrix_set(M, i, j, gsl_matrix_get(M, i, j) - tmp);
  }

  /* multiply by -1/2 */
  for (i = 0; i < nb_indiv; i++){
    for (j = 0; j < nb_indiv; j++){
      gsl_matrix_set(M, i, j, -0.5 * gsl_matrix_get(M, i, j));
    }
  }

  /* compute eigenvectors and eigenvalues */
  eval = gsl_vector_alloc(nb_indiv);
  evec = gsl_matrix_alloc(nb_indiv, nb_indiv);
  w = gsl_eigen_symmv_alloc(nb_indiv);
  gsl_eigen_symmv(M, eval, evec, w);
  gsl_eigen_symmv_free(w);
  gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_DESC);

  /* final step and put everything in double** to simplify next steps */
  mds = allocation_double_matrix(nb_indiv, NB_DIM_MDS);
  for (i = 0; i < nb_indiv; i++){
    for (j = 0; j < NB_DIM_MDS; j++){
      mds[i][j] = gsl_matrix_get(evec, i, j) * sqrt(gsl_vector_get(eval, j));
    }
  }
  
  gsl_matrix_free(M);
  gsl_matrix_free(evec);
  gsl_vector_free(eval);
  return mds;
}


/* computes the ASD matrix using multithreading over SNPs */
/* gives the same results as https://github.com/szpiech/asd */
double **compute_asd_matrix_local(param *params, arg *args, char **genotypes){
  double **asd = NULL;
  unsigned i, j, k, increment;
  pthread_t *threads = NULL;
  thread_datast_asd *thread_data = NULL;
  int dummy;
  unsigned nb_indiv = args->sample_size_s1 + args->sample_size_s2 + args->sample_size_adm;
  fprintf(stderr, "Computing ASD matrix...\n");
  
  threads = allocation_pthread_t_vector(args->nb_thread);
  thread_data = allocation_thread_datast_asd_vector(args->nb_thread);

  increment = args->nb_snp / args->nb_thread;
  for (i = 0; i < args->nb_thread; i++){
    thread_data[i].tmp_asd = allocation_double_matrix(nb_indiv, nb_indiv);
    thread_data[i].genotypes = genotypes;
    thread_data[i].idx_min = increment * i;
    thread_data[i].idx_max = increment * (i+1) - 1;
    if (i == args->nb_thread - 1)
      thread_data[i].idx_max = args->nb_snp - 1;
    thread_data[i].nb_indiv = nb_indiv;
  }
  for (k = 0; k < args->nb_thread; k++){
    dummy = pthread_create(&threads[k], NULL, compute_asd_snp_subset, (void*) &thread_data[k]);
    if (dummy != 0)
      exit_on_error(NULL, PTHREAD_ERROR);
  }    
  for (k = 0; k < args->nb_thread; k++)
    pthread_join(threads[k], NULL);

  asd = allocation_double_matrix(nb_indiv, nb_indiv);
  for (i = 0; i < nb_indiv; i++){
    asd[i][i] = 0.0;
    for (j = i+1; j < nb_indiv; j++){
      asd[i][j] = 0.0;
      for (k = 0; k < args->nb_thread; k++)
	asd[i][j] += thread_data[k].tmp_asd[i][j];

      asd[i][j] /= args->nb_snp;
      asd[j][i] = asd[i][j];
    }
  }
  for (k = 0; k < args->nb_thread; k++)
    free_double_matrix(thread_data[k].tmp_asd, nb_indiv);

  free(thread_data);
  free(threads);
  
  return asd;
}

/* computes distance on a subset of SNPs */
void *compute_asd_snp_subset(void *ptr){
  thread_datast_asd *data = (thread_datast_asd*) ptr;
  unsigned i, j, k;
  for (i = 0; i < data->nb_indiv; i++){
    data->tmp_asd[i][i] = 0.0;
    for (j = i+1; j < data->nb_indiv; j++){
      data->tmp_asd[i][j] = 0.0;
      for (k = data->idx_min; k <= data->idx_max; k++){
	if (data->genotypes[i][k] == data->genotypes[j][k]){	/* all alleles shared -> dist=0 */
	  continue;
	} else {
	  if (data->genotypes[i][k] == 1 || data->genotypes[j][k] == 1){
	    data->tmp_asd[i][j] += 0.5;	/* 1 vs * -> .5 allele shared */
	  } else {
	    data->tmp_asd[i][j] += 1.0;
	  }
	}
      }
    }
  }
  
  return NULL;
}

/* changes the format of genotypes for simpler sumstats computations */
/* multithreaded over SNPs */
char **simplify_genotypes(param *params, arg *args, unsigned **wanted_samples){
  char **genotypes = NULL;
  pthread_t *threads = NULL;
  thread_datast_fill_genotypes *thread_datas = NULL;
  unsigned nb_indiv = args->sample_size_s1 + args->sample_size_s2 + args->sample_size_adm;
  unsigned i, increment, idx_deb, idx_fin;
  int dummy;
  
  fprintf(stderr, "Simplifying genotypes matrix...\n");
  
  threads = allocation_pthread_t_vector(args->nb_thread);
  thread_datas = allocation_thread_datast_fill_genotypes_vector(args->nb_thread);
  genotypes = allocation_char_matrix(nb_indiv, args->nb_snp);
  increment = args->nb_snp / args->nb_thread;
  for (i = 0; i < args->nb_thread; i++){
    idx_deb = increment * i;
    idx_fin = increment * (i+1) - 1;
    if (i == args->nb_thread - 1)
      idx_fin = args->nb_snp - 1;
    thread_datas[i].idx_min = idx_deb;
    thread_datas[i].idx_max = idx_fin;
    thread_datas[i].params = params;
    thread_datas[i].args = args;
    thread_datas[i].wanted_samples = wanted_samples;
    thread_datas[i].genotypes = genotypes;
    dummy = pthread_create(&threads[i], NULL, fill_genotype_matrix, (void*) &thread_datas[i]);
    if (dummy != 0)
      exit_on_error(NULL, PTHREAD_ERROR);
  }
  for (i = 0; i < args->nb_thread; i++){
    pthread_join(threads[i], NULL);
  }

  free(threads);
  free(thread_datas);
  
  return genotypes;
}

/* fills the genotypes matrix in a multithreaded way */
void *fill_genotype_matrix(void *ptr){
  thread_datast_fill_genotypes *data = (thread_datast_fill_genotypes*) ptr;
  unsigned i, j, done;
  char geno1, geno2, new_geno;
  for (j = data->idx_min; j <= data->idx_max; j++){
    done = 0;
    for (i = 0; i < data->args->sample_size_adm; i++){
      if (data->wanted_samples[done][0] != FLAG_ADM)
	printf("probleme...\n");
      geno1 = data->params->genos_adm_new[data->wanted_samples[done][1]][j];
      geno2 = data->params->genos_adm_new[data->wanted_samples[done][2]][j];
      new_geno = recode_geno(geno1, geno2);
      data->genotypes[i][j] = new_geno;
      done++;
    }
    for (i = 0; i < data->args->sample_size_s1; i++){
      if (data->wanted_samples[done][0] != FLAG_S1)
	printf("probleme...\n");
      geno1 = data->params->genos_s1[data->wanted_samples[done][1]][data->params->wanted_snps[j]];
      geno2 = data->params->genos_s1[data->wanted_samples[done][2]][data->params->wanted_snps[j]];
      new_geno = recode_geno(geno1, geno2);
      data->genotypes[done][j] = new_geno;
      done++;
    }
    for (i = 0; i < data->args->sample_size_s2; i++){
      if (data->wanted_samples[done][0] != FLAG_S2)
	printf("probleme...\n");
      geno1 = data->params->genos_s2[data->wanted_samples[done][1]][data->params->wanted_snps[j]];
      geno2 = data->params->genos_s2[data->wanted_samples[done][2]][data->params->wanted_snps[j]];
      new_geno = recode_geno(geno1, geno2);
      data->genotypes[done][j] = new_geno;
      done++;
    }
  }
  return NULL;
}

char recode_geno(char geno1, char geno2){
  if (geno1 != geno2)
    return 1;
  if (geno1 == '0')
    return 0;
  return 2;
}

thread_datast_fill_genotypes *allocation_thread_datast_fill_genotypes_vector(unsigned nb_thread){
  thread_datast_fill_genotypes *v = NULL;
  v = (thread_datast_fill_genotypes*) malloc(nb_thread * sizeof(thread_datast_fill_genotypes));
  check_alloc(v);
  return v;
}

thread_datast_asd *allocation_thread_datast_asd_vector(unsigned nb_thread){
  thread_datast_asd *v = NULL;
  v = (thread_datast_asd*) malloc(nb_thread * sizeof(thread_datast_asd));
  check_alloc(v);
  return v;
}
