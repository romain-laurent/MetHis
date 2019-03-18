#include "misc.h"
#include "arguments_parser.h"
#include <getopt.h>
#include <errno.h>
#include <ctype.h>
#include <string.h>

arg *parse_arguments(int argc, char **argv){
  arg *args = NULL;
  int c, option_index, dummy;
  char *dummy2, *dummy3, *dummy4;
  dummy2 = allocation_char_vector(SMALL_BUFF_SIZE);
  dummy3 = allocation_char_vector(SMALL_BUFF_SIZE);
  dummy4 = allocation_char_vector(SMALL_BUFF_SIZE);

  args = allocation_arg();
  
  if (argc == 1){
    args->help = 1;
  }
  
  struct option long_options[] = {
    /* flag */
    {"help", no_argument, &(args->help), 1},
    {"save_data", no_argument, &(args->save_data), 1},
    {"force_rewrite", no_argument, &(args->force_rewrite), 1},
    {"preview", no_argument, &(args->draw_preview), 1},
    /* require arguments */
    {"nb_snp", required_argument, 0, 'a'},
    {"nb_simul", required_argument, 0, 'c'},
    {"nb_thread", required_argument, 0, 'd'},
    {"sampling", required_argument, 0, 'e'},
    {"input_path", required_argument, 0, 'f'},
    {"prefix", required_argument, 0, 'g'},
    {"plots", required_argument, 0, 'l'},
    {0, 0, 0, 0}
  };
  
  while (1 && args->help != 1){
    option_index = 0;
    c = getopt_long(argc, argv, "a:c:d:e:f:g:l:", long_options, &option_index);
    
    /* if end of options */
    if (c == -1)
      break;
    
    /* unrecognized option */
    if (c == '?' || c == ':'){
      fprintf(stderr, "Unrecognized option. Exiting\n");
      exit_on_error(NULL, NO_CODE);
    }
    
    /* else, we deal with options */
    switch (c){
    case 0 :
      break;
    case 'a' :
      args->nb_snp = get_unsigned_number(optarg, "nb_snp");
      break;
    case 'c' :
      dummy = sscanf(optarg, "%[^-]-%s", dummy2, dummy3);
      if (dummy == 2){
	args->idx_simul_deb = get_unsigned_number(dummy2, "nb_simul");
	args->idx_simul_fin = get_unsigned_number(dummy3, "nb_simul");
      }
      else{
	dummy = sscanf(optarg, "%s", dummy2);
	args->idx_simul_deb = 1;
	args->idx_simul_fin = get_unsigned_number(dummy2, "nb_simul");
      }
      break;
    case 'd' :
      args->nb_thread = get_unsigned_number(optarg, "nb_thread");
      break;
    case 'e' :
      dummy = sscanf(optarg, "%[^/]/%[^/]/%s", dummy2, dummy3, dummy4);
      args->sample_size_s1 = get_unsigned_number(dummy2, "sampling");
      args->sample_size_adm = get_unsigned_number(dummy3, "sampling");
      args->sample_size_s2 = get_unsigned_number(dummy4, "sampling");
      break;
    case 'f':
      dummy = sscanf(optarg, "%s", args->input_path);
      break;
    case 'g' :
      dummy = sscanf(optarg, "%s", args->prefix);
      break;
    case 'l':
      args->draw_plots = parse_plots_argument(optarg);
      break;
    }
  }

  free(dummy2);
  free(dummy3);
  free(dummy4);
  
  if (args->help)
    print_usage(argv[0]);
  check_args(args);
  return args;
}

void check_args(arg *args){
  /* first we check easy arguments */
  if (args->nb_snp == 0){
    fprintf(stderr, "Error: nb_snp cannot be 0. Exiting\n");
    exit_on_error(NULL, NO_CODE);
  }
  if (args->nb_thread == 0){
    fprintf(stderr, "Error: nb_thread cannot be 0. Exiting.\n");
    exit_on_error(NULL, NO_CODE);
  }
  if (args->idx_simul_deb > args->idx_simul_fin){
    fprintf(stderr, "Error: invalid argument for nb_simul. Exiting.\n");
    exit_on_error(NULL, NO_CODE);
  }
}

unsigned parse_plots_argument(char *option){
  char *tmp = NULL;
  unsigned val = 0;
  tmp = allocation_char_vector(SMALL_BUFF_SIZE);
  tmp = strncpy(tmp, option, SMALL_BUFF_SIZE);
  if (strncmp(tmp, "None", SMALL_BUFF_SIZE) == 0)
    val = PLOT_NONE;
  else if (strncmp(tmp, "All_gen", SMALL_BUFF_SIZE) == 0)
    val = PLOT_ALL;
  else if (strncmp(tmp, "Last_gen", SMALL_BUFF_SIZE) == 0)
    val = PLOT_LAST;
  else{
    fprintf(stderr, "Invalid argument for plots\n");
    exit_on_error(NULL, NO_CODE);
  }
  free(tmp);
  return val;
}

unsigned get_unsigned_number(char *str, char *option_name){
  unsigned val;
  char *endptr;
  if (! is_only_digits(str)){
    fprintf(stderr, "Invalid argument for %s\n", option_name);
    exit_on_error(NULL, NO_CODE);
  }
  errno = 0;
  val = (unsigned) strtoul(str, &endptr, 10);
  if (errno != 0 || endptr == str){
    fprintf(stderr, "Invalid argument for %s\n", option_name);
    exit_on_error(NULL, NO_CODE);
  }
  return val;
}

unsigned is_only_digits(char *str){
  unsigned i;
  for (i = 0; str[i] != '\0'; i++)
    if (! isdigit(str[i]))
      return 0;
  return 1;
}

	
