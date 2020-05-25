# your choice of compiler
CC = mpicc

# Add your choice of flags
CFLAGS = -O3 -Wall -Wextra -g -fopenmp
LDLIBS = -lm
MATRIX = cfd1.mtx

all : cg_mpi checker

cg_mpi : cg_mpi.o mmio.o
	$(CC) $(CFLAGS) -o cg_mpi cg_mpi.o mmio.o $(LDLIBS)

checker : checker.o mmio.o
	gcc $(CFLAGS) -o checker checker.o mmio.o $(LDLIBS)

mmio.o : mmio.c mmio.h
checker.o : checker.c mmio.h
cg_mpi.o : cg_mpi.c mmio.h

exec : cg_mpi
	 mpiexec -n 3 sh -c "zcat ${MATRIX}.gz | ./$< --seed 23 --solution x.txt"

execCheck : checker
	sh -c "zcat ${MATRIX}.gz | ./$< --seed 23 --solution x.txt"

.PHONY: clean
clean :
	rm -rf *.o cg_mpi checker
