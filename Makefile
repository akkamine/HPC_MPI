# your choice of compiler
CC = gcc

# Add your choice of flags
CFLAGS = -O3 -Wall -Wextra -pg -g -fopenmp
LDLIBS = -lm

all : cg

cg : cg.o mmio.o
		$(CC) $(CFLAGS) -o cg cg.o mmio.o $(LDLIBS)

mmio.o : mmio.c mmio.h
cg.o : cg.c mmio.h

.PHONY: clean
clean :
	rm -rf cg *.o
