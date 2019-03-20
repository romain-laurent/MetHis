
#ifndef SET_PARAMS_H
#define SET_PARAMS_H

unsigned check_param_file_nb_field(FILE *f);
void read_param_file(param *params, arg *args, unsigned idx_simul);
void seed_RNGs_from_file(param *params, unsigned nb_rng);
unsigned get_seed_from_file(void);
void write_next_seed(unsigned seed);
void create_simul_log(arg *args, param *params, unsigned idx_simul);
void reset_param(param *params, arg *args);
void set_params_one_simul(param *params, arg *args, unsigned idx_simul);
void reseed_rngs(param *params, unsigned nb_rng, char *mess);

#endif
