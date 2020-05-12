CC = gcc
CC_MPI = mpicc

# Add your choice of flags
CFLAGS = -O3 -w -g -fopenmp -ftree-vectorize -mavx2 -mfma
MPIFLAGS = -O3 -w -g -ftree-vectorize -mavx2 -mfma
LDLIBS = -lm
TARGET = cg_mpi
MATRIX = tmt_sym.mtx


all : cg_mpi

exec_mpi : cg_mpi
	mpiexec --n 2 --display-map ./$< --matrix ${MATRIX} > /dev/null

#object files
mmio.o : mmio.c mmio.h

cg_mpi.o : cg_mpi.c mmio.h
	${CC_MPI} ${MPIFLAGS} -c -o $@ $< ${LDLIBS}

%.o : %.c mmio.h
	${CC} ${CFLAGS} -c -o $@ $< ${LDLIBS}

#executable
cg_mpi : cg_mpi.o mmio.h
	${CC_MPI} ${CFLAGS} -o $@ $^ ${LDLIBS}

.PHONY: clean
clean :
	rm -rf *.o ${TARGET}
