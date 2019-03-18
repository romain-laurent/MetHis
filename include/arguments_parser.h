
#ifndef ARGUMENTS_PARSER_H
#define ARGUMENTS_PARSER_H

arg *parse_arguments(int argc, char **argv);
void parse_Ne_argument(char *option, arg *args);
unsigned is_only_digits(char *v);
void parse_contrib_argument(char *option, arg *args, char id);
unsigned get_unsigned_number(char *str, char *option_name);
unsigned parse_model(char *str, char *option_name);
unsigned parse_plots_argument(char *option);
double get_double_number(char *str, char id);
void parse_contrib_trend(arg *args, char *option, char id);
void parse_contrib_x(arg *args, char *option, char id, char id2);
void check_args(arg *args);


#endif
