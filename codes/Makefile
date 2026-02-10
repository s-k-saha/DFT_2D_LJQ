CC      = gcc
CFLAGS = -Wall -Wextra -std=c11 -O0 -g
LDLIBS  = -lm -lfftw3 

MAIN    = codeLJQ.c
UTIL    = utility.c
OBJ     = $(MAIN:.c=.o) $(UTIL:.c=.o)

EXEC    = ../exec/$(basename $(notdir $(MAIN))).out

.PHONY: all clean

all: $(EXEC)

$(EXEC): $(OBJ)
	@mkdir -p ../exec
	$(CC) -o $@ $(OBJ) $(LDLIBS)

%.o: %.c utility.h
	$(CC) -c $<

clean:
	rm -f $(OBJ) $(EXEC)

