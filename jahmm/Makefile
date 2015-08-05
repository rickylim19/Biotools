SRC_DIR= src
INC_DIR= src
LIB_DIR= lib
OBJECT_FILES= jahmm.o zinb.o hmm.o utils.o xxhash.o
SOURCE_FILES= main.c samread.c

OBJECTS= $(addprefix $(SRC_DIR)/,$(OBJECT_FILES))
SOURCES= $(addprefix $(SRC_DIR)/,$(SOURCE_FILES))
INCLUDES= $(addprefix -I, $(INC_DIR))

P= jahmm
CC= gcc
#CFLAGS= -std=gnu99 -g -Wall -O0
CFLAGS= -std=gnu99 -g -Wall -O3
LDLIBS= -lm -Llib -Wl,-rpath=lib -lhts -lsam

all: $(P)

$(P): $(OBJECTS) $(SOURCES)
	$(CC) $(CFLAGS) $(OBJECTS) $(SOURCES) $(LDLIBS) -o $@

$(SRC_DIR)/%.o: $(SRC_DIR)/%.c $(SRC_DIR)/%.h
	        $(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(P)
