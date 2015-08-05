SRC_DIR= ../../src
INC_DIR= ../../src

INCLUDES= $(addprefix -I, $(INC_DIR))
OBJECTS= jahmm.o zinb.o hmm.o utils.o xxhash.o

CC= gcc
CFLAGS= -std=gnu99 -fPIC -g -Wall -O3

all: Rjahmm.so

%.o: $(SRC_DIR)/%.c $(SRC_DIR)/%.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

Rjahmm.so: $(OBJECTS)
	R CMD SHLIB $(INCLUDES) Rjahmm.c $(OBJECTS)

clean:
	rm -f $(OBJECTS) Rjahmm.o Rjahmm.so
