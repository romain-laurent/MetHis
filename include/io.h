
#ifndef IO_H
#define IO_H

#define t_or_f(x) ((x) == 1 ? "True" : "False")

void print_number_or_NA(FILE *f, double N);
void print_sumstats(sumstats_datast *sumstats, arg *args, param *params);
void write_arguments(arg *args, FILE *log);
void get_time(char *str);
void  read_input_genotypes(param *params, arg *args);
unsigned count_nb_chrom(FILE *f, char *buff, char id);
unsigned count_nb_snp(char *str);
unsigned skip_next_space(char *str, unsigned idx);
void print_params_ranges(param *params, arg *args);
void write_vcf_file(param *params, arg *args, unsigned **wanted_samples, unsigned idx_gen);


#endif
