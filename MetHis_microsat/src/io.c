#include "misc.h"
#include "io.h"
#include <time.h>
#include <string.h>
#include <ctype.h>

/* Input/ouput fonctions */

/* creates admixtools inputs using bash and vcftools */
void create_admixtools_inputs(param *params, arg *args, unsigned idx_gen){
  char *tmp = NULL, *command = NULL;
  FILE *f = NULL;
  int dummy;
  command = allocation_char_vector(LARGE_BUFF_SIZE);
  tmp = allocation_char_vector(MEDIUM_BUFF_SIZE);
  sprintf(tmp, "%s/simu_%u/simu_%u_g%u", args->prefix, params->current_simul, params->current_simul, idx_gen);
  /* create PED/MAP files */
  sprintf(command, "%s --vcf %s.vcf --plink --out %s 2> /dev/null", VCFTOOLS_PATH, tmp, tmp);
  dummy = system(command);
  if (dummy != 0){
    fprintf(stderr, "It seems something went wrong during generation of ADMIXTOOLS input\n");
    fprintf(params->simul_log, "It seems something went wrong during generation of ADMIXTOOLS input\n");
  }
  
  /* create ind file */
  sprintf(tmp, "%s/simu_%u", args->prefix, params->current_simul);
  sprintf(command, "for type in adm s1 s2; do sed \"s/$/ M $type/\" %s/indivs_$type.txt >> %s/indivs.ind; done", tmp, tmp);
  dummy = system(command);
  if (dummy != 0){
    fprintf(stderr, "It seems something went wrong during generation of ADMIXTOOLS input\n");
    fprintf(params->simul_log, "It seems something went wrong during generation of ADMIXTOOLS input\n");
  }

  /* create list file */
  sprintf(command, "%s/list_qp3pop.txt", tmp);
  f = safe_open(command, "w");
  fprintf(f, "s1 s2 adm\n");
  fclose(f);

  /* create param file */
  sprintf(command, "%s/qp3pop.param", tmp);
  f = safe_open(command, "w");
  fprintf(f, "indivname: %s/indivs.ind\n", tmp);
  fprintf(f, "popfilename: %s/list_qp3pop.txt\n", tmp);
  sprintf(tmp, "%s/simu_%u/simu_%u_g%u", args->prefix, params->current_simul, params->current_simul, idx_gen);
  fprintf(f, "snpname: %s.map\n", tmp);
  fprintf(f, "genotypename: %s.ped\n", tmp);
  fclose(f);
  
  free(command);
  free(tmp);
}


/* writes a vcf file and all we need to use vcftools */
void write_vcf_file(param *params, arg *args, unsigned **wanted_samples, unsigned idx_gen){
  FILE *f = NULL;
  char buff[LARGE_BUFF_SIZE];
  unsigned i, j, nb_done;
  /* first we create the files for each pop */
  sprintf(buff, "./%s/simu_%u/indivs_s1.txt", args->prefix, params->current_simul);
  f = safe_open(buff, "w");
  for (i = 0; i < args->sample_size_s1; i++)
    fprintf(f, "s1_%u\n", i);
  fclose(f);
  sprintf(buff, "./%s/simu_%u/indivs_s2.txt", args->prefix, params->current_simul);
  f = safe_open(buff, "w");
  for (i = 0; i < args->sample_size_s2; i++)
    fprintf(f, "s2_%u\n", i);
  fclose(f);
  sprintf(buff, "./%s/simu_%u/indivs_adm.txt", args->prefix, params->current_simul);
  f = safe_open(buff, "w");
  for (i = 0; i < args->sample_size_adm; i++)
    fprintf(f, "adm_%u\n", i);
  fclose(f);
  /* then we write the vcf */
  sprintf(buff, "./%s/simu_%u/simu_%u_g%u.vcf", args->prefix, params->current_simul, params->current_simul, idx_gen);
  f = safe_open(buff, "w");
  fprintf(f, "##fileformat=VCFv4.1\n");
  /* header */
  fprintf(f, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
  for (i = 0; i < args->sample_size_adm; i++)
    fprintf(f, "\tadm_%u", i);
  for (i = 0; i < args->sample_size_s1; i++)
    fprintf(f, "\ts1_%u", i);
  for (i = 0; i < args->sample_size_s2; i++)
    fprintf(f, "\ts2_%u", i);
  fprintf(f, "\n");
  /* genotypes */
  for (i = 0; i < args->nb_snp; i++){
    fprintf(f, "1\t%d\t.\tA\tCN0", (i+1)*2000);
    {
      unsigned tmp_i;
      for (tmp_i=1; tmp_i < 100; tmp_i++)
	fprintf(f, ",CN%u", tmp_i);
    }
    fprintf(f, "\t1\tPASS\tD\tGT");
    nb_done = 0;
    for (j = 0; j < args->sample_size_adm; j++){
      if (wanted_samples[nb_done][0] != FLAG_ADM)
	printf("probleme...\n");
      fprintf(f, "\t%d/%d", params->genos_adm_new[wanted_samples[nb_done][1]][i], params->genos_adm_new[wanted_samples[nb_done][2]][i]);
      nb_done++;
    }
    for (j = 0; j < args->sample_size_s1; j++){
      if (wanted_samples[nb_done][0] != FLAG_S1)
	printf("probleme...\n");
      fprintf(f, "\t%d/%d", params->genos_s1[wanted_samples[nb_done][1]][params->wanted_snps[i]], params->genos_s1[wanted_samples[nb_done][2]][params->wanted_snps[i]]);
      nb_done++;
    }
    for (j = 0; j < args->sample_size_s2; j++){
      if (wanted_samples[nb_done][0] != FLAG_S2)
	printf("probleme...\n");
      fprintf(f, "\t%d/%d", params->genos_s2[wanted_samples[nb_done][1]][params->wanted_snps[i]], params->genos_s2[wanted_samples[nb_done][2]][params->wanted_snps[i]]);
      nb_done++;
    }
    fprintf(f, "\n");
  }
  fclose(f);
}


/* write a tped (and a tfam) file */
void write_tped_file(param *params, arg *args, unsigned **wanted_samples, unsigned idx_gen){
  FILE *f = NULL;
  char buff[MEDIUM_BUFF_SIZE];
  unsigned i, j, nb_done;
  /* first we create the tfam file */
  sprintf(buff, "./%s/simu_%u/simu_%u_g%u.tfam", args->prefix, params->current_simul, params->current_simul, idx_gen);
  f = safe_open(buff, "w");
  for (i = 0; i < args->sample_size_adm; i++)
    fprintf(f, "adm_%u\tadm_%u\t0\t0\t-9\n", i, i);
  for (i = 0; i < args->sample_size_s1; i++)
    fprintf(f, "s1_%u\ts1_%u\t0\t0\t-9\n", i, i);
  for (i = 0; i < args->sample_size_s2; i++)
    fprintf(f, "s2_%u\ts2_%u\t0\t0\t-9\n", i, i);
  fclose(f);
  /* then the tped */
  sprintf(buff, "./%s/simu_%u/simu_%u_g%u.tped", args->prefix, params->current_simul, params->current_simul, idx_gen);
  f = safe_open(buff, "w");
  for (i = 0; i < args->nb_snp; i++){
    fprintf(f, "%u\tsnp%u\t0\t1", i+1, i+1);
    nb_done = 0;
    for (j = 0; j < args->sample_size_adm; j++){
      if (wanted_samples[nb_done][0] != FLAG_ADM)
	printf("probleme...\n");
      fprintf(f, "\t%d\t%d", params->genos_adm_new[wanted_samples[nb_done][1]][i]+1, params->genos_adm_new[wanted_samples[nb_done][2]][i]+1);
      nb_done++;
    }
    for (j = 0; j < args->sample_size_s1; j++){
      if (wanted_samples[nb_done][0] != FLAG_S1)
	printf("probleme...\n");
      fprintf(f, "\t%d\t%d", params->genos_s1[wanted_samples[nb_done][1]][params->wanted_snps[i]]+1, params->genos_s1[wanted_samples[nb_done][2]][params->wanted_snps[i]]+1);
      nb_done++;
    }
    for (j = 0; j < args->sample_size_s2; j++){
      if (wanted_samples[nb_done][0] != FLAG_S2)
	printf("probleme...\n");
      fprintf(f, "\t%d\t%d", params->genos_s2[wanted_samples[nb_done][1]][params->wanted_snps[i]]+1, params->genos_s2[wanted_samples[nb_done][2]][params->wanted_snps[i]]+1);
      nb_done++;
    }
    fprintf(f, "\n");
  }
  fclose(f);
}



/* function to parse Arlequin file */
/* may not be super robust... */
/* assumes that fsc files will always be formatted as the example file provided by Cesar */
void  read_input_genotypes(param *params, arg *args){
  FILE *f_in = safe_open(args->input_path, "r");
  char *buff = allocation_char_vector(LARGE_BUFF_SIZE);
  unsigned pop_started = 0;
  char wanted_pop = '1';
  unsigned idx_buff, idx_snp, idx_chrom = 0;
  int geno, dummy;
  /* first we count the number of chromosome for each source population */
  params->nb_chrom_s1 = count_nb_chrom(f_in, buff, '1');
  params->nb_chrom_s2 = count_nb_chrom(f_in, buff, '2');
  /* check that we will have enough chromosomes for the sample size required */
  if (params->nb_chrom_s1 < 2 * args->sample_size_s1){
    sprintf(buff, "Sampling size error for source population s1: %u chromosomes required but only %u in input file. Exiting\n", 2 * args->sample_size_s1, params->nb_chrom_s1);
    exit_on_error(buff, SSIZE_ERROR);
  }
  if (params->nb_chrom_s2 < 2 * args->sample_size_s2){
    sprintf(buff, "Sampling size error for source population s2: %u chromosomes required but only %u in input file. Exiting\n", 2 * args->sample_size_s2, params->nb_chrom_s2);
    exit_on_error(buff, SSIZE_ERROR);
  }
  /* then we get the genotypes */
  while (fgets(buff, LARGE_BUFF_SIZE, f_in) != NULL){
    /* if we finished reading genotypes for one pop */
    if (pop_started && buff[0] == '}'){
      idx_chrom = 0;
      wanted_pop++;
      continue;
    }
    /* if the line does not contain genotypes */
    if (buff[0] != wanted_pop)
      continue;
    /* now there are some genotypes */
    /* if it's the first time we see genotypes, we count the number of SNPs and allocate memory to store them */
    if (params->genos_s1 == NULL){
      params->total_nb_snp = count_nb_snp(buff);
      params->genos_s1 = allocation_int_matrix(params->nb_chrom_s1, params->total_nb_snp);
      params->genos_s2 = allocation_int_matrix(params->nb_chrom_s2, params->total_nb_snp);
    }
    pop_started = 1;
    idx_buff = skip_next_space(buff, 0);
    idx_buff = skip_next_space(buff, idx_buff);
    idx_snp = 0;
    while (buff[idx_buff] != '\n'){
      if (buff[idx_buff] == '\0' || idx_snp >= params->total_nb_snp)
	exit_on_error(NULL, PARSING_ERROR);
      if ((wanted_pop == '1' && idx_chrom >= params->nb_chrom_s1) || (wanted_pop == '2' && idx_chrom >= params->nb_chrom_s2))
	exit_on_error(NULL, PARSING_ERROR);

      dummy = sscanf(buff + idx_buff, "%d", &geno);
      if (dummy != 1)
	exit_on_error(NULL, PARSING_ERROR);

      if (wanted_pop == '1')
	params->genos_s1[idx_chrom][idx_snp] = geno;
      else if (wanted_pop == '2')
	params->genos_s2[idx_chrom][idx_snp] = geno;
      idx_buff = skip_next_space(buff, idx_buff);
      idx_snp++;
      if (buff[idx_buff] == 0)
	break;
    }
    if (idx_snp != params->total_nb_snp)
      exit_on_error(NULL, PARSING_ERROR);
    idx_chrom++;
  }
  if (params->total_nb_snp < args->nb_snp)
    exit_on_error(NULL, NB_SNP_ERROR);
  free(buff);
  fclose(f_in);


}

unsigned count_nb_snp(char *str){
  unsigned i = 1;
  unsigned nb_snp = 0;
  /* counts number of spaces */
  while (str[i] != '\n'){
    if ( (!(isspace(str[i]))) && isspace(str[i-1]))
      nb_snp++;
    i++;
  }
  /* there is one less SNP */
  nb_snp--;
  return nb_snp;
}

unsigned skip_next_space(char *str, unsigned idx){
  unsigned space_found = 0;
  while (str[idx] != '\0'){
    if ((! space_found) && isspace(str[idx])){
      space_found = 1;
      idx++;
      continue;
    }
    if (space_found && (! isspace(str[idx])))
      break;
    idx++;
  }
  if (! space_found)
    exit_on_error(NULL, PARSING_ERROR);
  return idx;
}

unsigned count_nb_chrom(FILE *f, char *buff, char id){
  unsigned nb_chrom = 0;
  while(fgets(buff, LARGE_BUFF_SIZE, f) != NULL){
    if (buff[0] == id)
      nb_chrom++;
  }
  rewind(f);
  return nb_chrom;
}

/* writes arguments interpretation to both stderr and the general log */
void write_arguments(arg *args, FILE *log){
  char *line = NULL, *line2 = NULL;
  line = allocation_char_vector(MEDIUM_BUFF_SIZE);
  line2 = allocation_char_vector(MEDIUM_BUFF_SIZE);
  sprintf(line, "This is %s.\n", PROG_NAME);
  fprintf(stderr, "%s", line);
  fprintf(log, "%s", line);
  get_time(line);
  fprintf(stderr, "%s", line);
  fprintf(log, "%s", line);
  fprintf(stderr, "Command line interpreted as:\n");
  fprintf(log, "Command line interpreted as:\n");
  fprintf(stderr, "\tsave-data: %s\n", t_or_f(args->save_data));
  fprintf(log, "\tsave-data: %s\n", t_or_f(args->save_data));
  fprintf(stderr, "\tnb-snp: %u\n", args->nb_snp);
  fprintf(log, "\tnb-snp: %u\n", args->nb_snp);
  fprintf(stderr, "\tnb-simul: %u-%u\n", args->idx_simul_deb, args->idx_simul_fin);
  fprintf(log, "\tnb-simul: %u-%u\n", args->idx_simul_deb, args->idx_simul_fin);
  fprintf(stderr, "\tnb-thread: %u\n", args->nb_thread);
  fprintf(log, "\tnb-thread: %u\n", args->nb_thread);
  fprintf(stderr, "\tsampling: %u/%u/%u\n", args->sample_size_s1, args->sample_size_adm, args->sample_size_s2);
  fprintf(log, "\tsampling: %u/%u/%u\n", args->sample_size_s1, args->sample_size_adm, args->sample_size_s2);
  /* now options which are a mess... */
  fprintf(stderr, "\tinput-path: %s\n", args->input_path);
  fprintf(log, "\tinput-path: %s\n", args->input_path);
  free(line);
  free(line2);
}

/* stores the local time in a string in a human readable format */
void get_time(char *str){
  time_t t = time(NULL);
  struct tm tm = *localtime(&t);
  sprintf(str, "Local time is %02d-%02d-%d %02d:%02d:%02d\n", tm.tm_mday, tm.tm_mon + 1, tm.tm_year + 1900, tm.tm_hour, tm.tm_min, tm.tm_sec);
}
