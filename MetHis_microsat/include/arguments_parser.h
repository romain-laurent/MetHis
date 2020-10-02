
#ifndef ARGUMENTS_PARSER_H
#define ARGUMENTS_PARSER_H


arg *parse_arguments(int argc, char **argv);
unsigned is_only_digits(char *v);
unsigned get_unsigned_number(char *str, char *option_name);
void check_args(arg *args);


#endif
