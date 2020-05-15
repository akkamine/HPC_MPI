# your choice of compiler
CC = mpicc

# Add your choice of flags
CFLAGS = -O3 -Wall -Wextra -g
LDLIBS = -lm
MATRIX = cfd1.mtx

all : cg_mpi

cg_mpi : cg_mpi.o mmio.o
	$(CC) $(CFLAGS) -o cg_mpi cg_mpi.o mmio.o $(LDLIBS)

mmio.o : mmio.c mmio.h
cg_mpi.o : cg_mpi.c mmio.h

exec_mpi : cg_mpi
	mpiexec -n 2 ./$< --matrix ${MATRIX} --seed 23 > /dev/null

.PHONY: clean
clean :
	rm -rf *.o cg_mpi
