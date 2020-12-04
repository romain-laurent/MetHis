
TARGET = MetHis
SRC = ./src/
OBJ = ./obj/
INCLUDE = ./include/
CC = gcc
CFLAGS = -Wall -ansi -pedantic -g -I$(INCLUDE) --std=c99 -O2
LDFLAGS = $(shell pkg-config --libs gsl) -lpthread

# definition of object files
SOURCES = $(shell find $(SRC) -name "*.c")
OBJECTS = $(patsubst $(SRC)%.c, $(OBJ)%.o, $(SOURCES))

# linking command
$(TARGET) : $(OBJECTS)
	$(CC) -o $@ $^ $(LDFLAGS)

# dependencies
$(OBJ)main.o : $(INCLUDE)misc.h $(INCLUDE)arguments_parser.h $(INCLUDE)io.h $(INCLUDE)set_params.h
$(OBJ)arguments_parser.o : $(INCLUDE)misc.h $(INCLUDE)arguments_parser.h
$(OBJ)io.o : $(INCLUDE)misc.h $(INCLUDE)io.h $(INCLUDE)sumstats.h
$(OBJ)set_params.o : $(INCLUDE)misc.h $(INCLUDE)set_params.h $(INCLUDE)io.h
$(OBJ)simul.o : $(INCLUDE)misc.h $(INCLUDE)simul.h $(INCLUDE)sumstats.h $(INCLUDE)io.h
$(OBJ)sumstats.o : $(INCLUDE)misc.h $(INCLUDE)sumstats.h $(INCLUDE)io.h

# compiling command
$(OBJ)%.o : $(SRC)%.c
	$(CC) $(CFLAGS) -c $< -o $@

# special targets
.PHONY : all clean re

all : $(TARGET)

clean :
	rm -f $(SRC)*~ $(OBJ)* $(INCLUDE)*~ *~ $(TARGET)

re : clean all
