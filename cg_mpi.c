/*
 * Sequential implementation of the Conjugate Gradient Method.
 *
 * Authors : Lilia Ziane Khodja & Charles Bouillaguet
 *
 * v1.02 (2020-04-3)
 *
 * CHANGE LOG:
 *    v1.01 : fix a minor printing bug in load_mm (incorrect CSR matrix size)
 *    v1.02 : use https instead of http in "PRO-TIP"
 *
 * USAGE:
 * 	$ ./cg --matrix bcsstk13.mtx                # loading matrix from file
 *      $ ./cg --matrix bcsstk13.mtx > /dev/null    # ignoring solution
 *	$ ./cg < bcsstk13.mtx > /dev/null           # loading matrix from stdin
 *      $  zcat matrix.mtx.gz | ./cg                # loading gziped matrix from
 *      $ ./cg --matrix bcsstk13.mtx --seed 42      # changing right-hand side
 *      $ ./cg --no-check < bcsstk13.mtx            # no safety check
 *
 * PRO-TIP :
 *      # downloading and uncompressing the matrix on the fly
 *	$ curl --silent https://hpc.fil.cool/matrix/bcsstk13.mtx.gz | zcat | ./cg
 */
#include <stdio.h>
#include <stdlib.h>
#include <err.h>
#include <math.h>
#include <getopt.h>
#include <sys/time.h>
#include <mpi.h>
#include <stddef.h>
#include <string.h>

#include "mmio.h"

#define THRESHOLD 1e-8		// maximum tolerance threshold

struct csr_matrix_t {
	int n;			// dimension
	int nz;			// number of non-zero entries
	int *Ap;		// row pointers
	int *Aj;		// column indices
	double *Ax;		// actual coefficient
};
/*************************** Utility functions ********************************/

/* Seconds (wall-clock time) since an arbitrary point in the past */
double wtime()
{
	struct timeval ts;
	gettimeofday(&ts, NULL);
	return (double)ts.tv_sec + ts.tv_usec / 1e6;
}

/* Pseudo-random function to initialize b (rumors says it comes from the NSA) */
#define ROR(x, r) ((x >> r) | (x << (64 - r)))
#define ROL(x, r) ((x << r) | (x >> (64 - r)))
#define R(x, y, k) (x = ROR(x, 8), x += y, x ^= k, y = ROL(y, 3), y ^= x)
double PRF(int i, unsigned long long seed)
{
	unsigned long long y = i, x = 0xBaadCafe, b = 0xDeadBeef, a = seed;
	R(x, y, b);
	for (int i = 0; i < 31; i++) {
		R(a, b, i);
		R(x, y, b);
	}
	x += i;
	union { double d; unsigned long long l;	} res;
	res.l = ((x << 2) >> 2) | (((1 << 10) - 1ll) << 52);
	return 2 * (res.d - 1.5);
}

/*************************** Matrix IO ****************************************/

/* Load MatrixMarket sparse symetric matrix from the file descriptor f */
struct csr_matrix_t *load_mm(FILE * f)
{
	MM_typecode matcode;
	int n, m, nnz;

	/* -------- STEP 1 : load the matrix in COOrdinate format */
	double start = wtime();

	/* read the header, check format */
	if (mm_read_banner(f, &matcode) != 0)
		errx(1, "Could not process Matrix Market banner.\n");
	if (!mm_is_matrix(matcode) || !mm_is_sparse(matcode))
		errx(1, "Matrix Market type: [%s] not supported (only sparse matrices are OK)", mm_typecode_to_str(matcode));
	if (!mm_is_symmetric(matcode) || !mm_is_real(matcode))
		errx(1, "Matrix type [%s] not supported (only real symmetric are OK)", mm_typecode_to_str(matcode));
	if (mm_read_mtx_crd_size(f, &n, &m, &nnz) != 0)
		errx(1, "Cannot read matrix size");
	fprintf(stderr, "[IO] Loading [%s] %d x %d with %d nz in triplet format\n", mm_typecode_to_str(matcode), n, n, nnz);
	fprintf(stderr, "     ---> for this, I will allocate %.1f MByte\n", 1e-6 * (40.0 * nnz + 8.0 * n));

	/* Allocate memory for the COOrdinate representation of the matrix (lower-triangle only) */
	int *Ti = malloc(nnz * sizeof(*Ti));
	int *Tj = malloc(nnz * sizeof(*Tj));
	double *Tx = malloc(nnz * sizeof(*Tx));
	if (Ti == NULL || Tj == NULL || Tx == NULL)
		err(1, "Cannot allocate (triplet) sparse matrix");

	/* Parse and load actual entries */
	for (int u = 0; u < nnz; u++) {
		int i, j;
		double x;
		if (3 != fscanf(f, "%d %d %lg\n", &i, &j, &x))
			errx(1, "parse error entry %d\n", u);
		Ti[u] = i - 1;	/* MatrixMarket is 1-based */
		Tj[u] = j - 1;
		/*
		 * Uncomment this to check input (but it slows reading)
		 * if (i < 1 || i > n || j < 1 || j > i)
		 *	errx(2, "invalid entry %d : %d %d\n", u, i, j);
		 */
		Tx[u] = x;
	}

	double stop = wtime();
	fprintf(stderr, "     ---> loaded in %.1fs\n", stop - start);

	/* -------- STEP 2: Convert to CSR (compressed sparse row) representation ----- */
	start = wtime();

	/* allocate CSR matrix */
	struct csr_matrix_t *A = malloc(sizeof(*A));
	if (A == NULL)
		err(1, "malloc failed");
	int *w = malloc((n + 1) * sizeof(*w));
	int *Ap = malloc((n + 1) * sizeof(*Ap));
	int *Aj = malloc(2 * nnz * sizeof(*Ap));
	double *Ax = malloc(2 * nnz * sizeof(*Ax));
	if (w == NULL || Ap == NULL || Aj == NULL || Ax == NULL)
		err(1, "Cannot allocate (CSR) sparse matrix");

	/* the following is essentially a bucket sort */

	/* Count the number of entries in each row */
	for (int i = 0; i < n; i++)
		w[i] = 0;
	for (int u = 0; u < nnz; u++) {
		int i = Ti[u];
		int j = Tj[u];
		w[i]++;
		if (i != j)	/* the file contains only the lower triangular part */
			w[j]++;
	}

	/* Compute row pointers (prefix-sum) */
	int sum = 0;
	for (int i = 0; i < n; i++) {
		Ap[i] = sum;
		sum += w[i];
		w[i] = Ap[i];
	}
	Ap[n] = sum;

	/* Dispatch entries in the right rows */
	for (int u = 0; u < nnz; u++) {
		int i = Ti[u];
		int j = Tj[u];
		double x = Tx[u];
		Aj[w[i]] = j;
		Ax[w[i]] = x;
		w[i]++;
		if (i != j) {	/* off-diagonal entries are duplicated */
			Aj[w[j]] = i;
			Ax[w[j]] = x;
			w[j]++;
		}
	}

	/* release COOrdinate representation */
	free(w);
	free(Ti);
	free(Tj);
	free(Tx);
	stop = wtime();
	fprintf(stderr, "     ---> converted to CSR format in %.1fs\n", stop - start);
	fprintf(stderr, "     ---> CSR matrix size = %.1fMbyte\n", 1e-6 * (24. * nnz + 4. * n));

	A->n = n;
	A->nz = sum;
	A->Ap = Ap;
	A->Aj = Aj;
	A->Ax = Ax;
	return A;
}

/*************************** Matrix accessors *********************************/

/* Copy the diagonal of A into the vector d. */
void extract_diagonal(const struct csr_matrix_t *A, double *d)
{
	int n = A->n;
	int *Ap = A->Ap;
	int *Aj = A->Aj;
	double *Ax = A->Ax;
	for (int i = 0; i < n; i++) {
		d[i] = 0.0;
		for (int u = Ap[i]; u < Ap[i + 1]; u++)
			if (i == Aj[u])
				d[i] += Ax[u];
	}
}

/* Matrix-vector product (with A in CSR format) : y = Ax */
void sp_gemv(const struct csr_matrix_t *A, const double *x, double *y)
{
	int n = A->n;
	int *Ap = A->Ap;
	int *Aj = A->Aj;
	double *Ax = A->Ax;
	for (int i = 0; i < n; i++) {
		y[i] = 0;
		for (int u = Ap[i]; u < Ap[i + 1]; u++) {
			int j = Aj[u];
			double A_ij = Ax[u];
			y[i] += A_ij * x[j];
		}
	}
}

/*************************** Vector operations ********************************/

/* dot product */
double dot(const int n, const double *x, const double *y)
{
	double sum = 0.0;
	for (int i = 0; i < n; i++)
		sum += x[i] * y[i];
	return sum;
}

/* euclidean norm (a.k.a 2-norm) */
double norm(const int n, const double *x)
{
	double sum = 0.0;
	for (int i = 0; i < n; i++)
		sum += x[i] * x[i];
	return sqrt(sum);
}

/*********************** conjugate gradient algorithm *************************/

/* Solve Ax == b (the solution is written in x). Scratch must be preallocated of size 6n */
void cg_solve(const struct csr_matrix_t *A, const double *b, double *x, const double epsilon, int my_rank, int nb_proc, double *scratch)
{
	fprintf(stderr, "rank == %d welcome in cg_solve n = %d, b[0] = %f, Ax[0]= %f \n", my_rank, A->n, b[0], A->Ax[0]);

	int n = A->n;
	int nz = A->nz;
	int nb_ligne = n / nb_proc;

	/*if(n % nb_proc == 0){
        //simplest case, same number of lines
        for(int i = 0; i < num_of_proc; i++){
            table1[i] = n/num_of_proc;
            table2[i] = (n/num_of_proc) * i;
        }
    }
    else{
        for(int i = 0; i < num_of_proc; i++){
            //+1 line for each of the first rem processes
            if(i < rem)
                table1[i] = (n/num_of_proc) +1;
            else
                table1[i] = n/num_of_proc;
            //calculate index of starting line accordingly
            if(i <= rem)
                table2[i] = ((n/num_of_proc) +1) * i;
            else
                table2[i] = (n/num_of_proc) * i + rem;
        }
    }
		*/

	if(my_rank == 0){

		fprintf(stderr, "[CG] Starting iterative solver\n");
		fprintf(stderr, "     ---> Working set : %.1fMbyte\n", 1e-6 * (12.0 * nz + 52.0 * n));
		fprintf(stderr, "     ---> Per iteration: %.2g FLOP in sp_gemv() and %.2g FLOP in the rest\n", 2. * nz, 12. * n);

		double *r = scratch;	        // residue
		double *z = scratch + n;	// preconditioned-residue
		double *p = scratch + 2 * n;	// search direction
		double *q = scratch + 3 * n;	// q == Ap
		double *d = scratch + 4 * n;	// diagonal entries of A (Jacobi preconditioning)

		/* Isolate diagonal */
		extract_diagonal(A, d);

		/*
		 * This function follows closely the pseudo-code given in the (english)
		 * Wikipedia page "Conjugate gradient method". This is the version with
		 * preconditionning.
		 */

		/* We use x == 0 --- this avoids the first matrix-vector product. */
		for (int i = 0; i < n; i++)
			x[i] = 0.0;
		for (int i = 0; i < n; i++)	// r <-- b - Ax == b
			r[i] = b[i];
		for (int i = 0; i < n; i++)	// z <-- M^(-1).r
			z[i] = r[i] / d[i];
		for (int i = 0; i < n; i++)	// p <-- z
			p[i] = z[i];

		double rz = dot(n, r, z);
		double start = wtime();
		double last_display = start;
		int iter = 0;
		while (norm(n, r) > epsilon) {
			/* loop invariant : rz = dot(r, z) */
			double old_rz = rz;
			sp_gemv(A, p, q);	/* q <-- A.p */
			double alpha = old_rz / dot(n, p, q);
			for (int i = 0; i < n; i++)	// x <-- x + alpha*p
				x[i] += alpha * p[i];
			for (int i = 0; i < n; i++)	// r <-- r - alpha*q
				r[i] -= alpha * q[i];
			for (int i = 0; i < n; i++)	// z <-- M^(-1).r
				z[i] = r[i] / d[i];
			rz = dot(n, r, z);	// restore invariant
			double beta = rz / old_rz;
			for (int i = 0; i < n; i++)	// p <-- z + beta*p
				p[i] = z[i] + beta * p[i];
			iter++;
			double t = wtime();
			if (t - last_display > 0.5) {
				/* verbosity */
				double rate = iter / (t - start);	// iterations per s.
				double GFLOPs = 1e-9 * rate * (2 * nz + 12 * n);
				fprintf(stderr, "\r     ---> error : %2.2e, iter : %d (%.1f it/s, %.2f GFLOPs)", norm(n, r), iter, rate, GFLOPs);
				fflush(stdout);
				last_display = t;
			}
		}
		fprintf(stderr, "\n     ---> Finished in %.1fs and %d iterations\n", wtime() - start, iter);
}
}

/******************************* main program *********************************/

/* options descriptor */
struct option longopts[6] = {
	{"seed", required_argument, NULL, 's'},
	{"rhs", required_argument, NULL, 'r'},
	{"matrix", required_argument, NULL, 'm'},
	{"solution", required_argument, NULL, 'o'},
	{"no-check", no_argument, NULL, 'c'},
	{NULL, 0, NULL, 0}
};

int main(int argc, char **argv)
{
	/* Parse command-line options */
	long long seed = 0;
	char *rhs_filename = NULL;
	char *matrix_filename = NULL;
	char *solution_filename = NULL;
	int safety_check = 1;
	char ch;
	while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
		switch (ch) {
		case 's':
			seed = atoll(optarg);
			break;
		case 'r':
			rhs_filename = optarg;
			break;
		case 'm':
			matrix_filename = optarg;
			break;
		case 'o':
			solution_filename = optarg;
			break;
		case 'c':
			safety_check = 0;
			break;
		default:
			errx(1, "Unknown option");
		}
	}

	int nb_proc;
	int my_rank;

	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	struct csr_matrix_t *A = malloc(sizeof(*A));
	struct csr_matrix_t *A2 = malloc(sizeof(*A2));
	double* x = NULL;
	double* b = NULL;
	double* y = NULL;
	double* scratch = NULL;

	int n = 0;
	int nz = 0;


	FILE *f_mat = stdin;
	if(my_rank == 0){
		/* Load the matrix */
		FILE *f_mat = stdin;
		if (matrix_filename) {
			f_mat = fopen(matrix_filename, "r");
			if (f_mat == NULL)
				err(1, "cannot matrix file %s", matrix_filename);
		}
		A = load_mm(f_mat);
		/* Allocate memory */
		n = A->n;
		nz = A->nz;
	}

	MPI_Bcast(&n , 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&nz , 1, MPI_INT, 0, MPI_COMM_WORLD);

	int *Ap = malloc(sizeof(int)*(n+1));
	int *Aj = malloc(sizeof(int)*nz);
	double *Ax = malloc(sizeof(double)*nz);

	if(my_rank ==0){
		memcpy(Ap, A->Ap, (n+1) * sizeof(int));
    memcpy(Aj, A->Aj, nz * sizeof(int));
    memcpy(Ax, A->Ax, nz * sizeof(double));
	}

	MPI_Bcast(Ap , n+1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(Aj , nz, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(Ax , nz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	fprintf(stderr, "rank = %d, n = %d, nz = %d, Ax[0]= %f\n", my_rank, n, nz, Ax[0]);
	MPI_Barrier(MPI_COMM_WORLD);

	A2->n = n;
	A2->nz = nz;
	A2->Ap = malloc(sizeof(int)*(n+1));
	A2->Aj = malloc(sizeof(int)*(nz));
	A2->Ax = malloc(sizeof(double)*(nz));

	memcpy(A2->Ap, Ap, (n+1) * sizeof(int));
	memcpy(A2->Aj, Aj, nz * sizeof(int));
	memcpy(A2->Ax, Ax, nz * sizeof(double));

	/*Creation du data_type csr_matrice pour broadcast la matrice à tout les processus*/
	/*MPI_Datatype mpi_csr_matrix;
	MPI_Datatype type[5] = { MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE };
	int blocklen[5] = { 1, 1, n+1, nz, nz };
	MPI_Aint disp[5];
  disp[0] = 0;
  disp[1] = disp[0] + blocklen[0] * sizeof(int);
  disp[2] = disp[1] + blocklen[1] * sizeof(int);
  disp[3] = disp[2] + blocklen[2] * sizeof(int);
  disp[4] = disp[3] + blocklen[3] * sizeof(int);

	MPI_Type_create_struct(5, blocklen, disp, type, &mpi_csr_matrix);
  MPI_Type_commit(&mpi_csr_matrix);

	if(my_rank ==0){
		tmp_A.n = A->n;
		tmp_A.nz = A->nz;
		memcpy(tmp_A.Ap, A->Ap, blocklen[2] * sizeof(int));
		memcpy(tmp_A.Aj, A->Aj, blocklen[3] * sizeof(int));
		memcpy(tmp_A.Ax, A->Ax, blocklen[4] * sizeof(double));
	}

	MPI_Bcast(&tmp_A , 1, mpi_csr_matrix, 0, MPI_COMM_WORLD);*/

	double *mem = malloc(7 * n * sizeof(double));
	if (mem == NULL)
		err(1, "cannot allocate dense vectors");
	x = mem;	/* solution vector */
	b = mem + n;	/* right-hand side */
	scratch = mem + 2 * n;	/* workspace for cg_solve() */

	if(my_rank == 0){
		/* Prepare right-hand size */
		if (rhs_filename) {	/* load from file */
			FILE *f_b = fopen(rhs_filename, "r");
			if (f_b == NULL)
				err(1, "cannot open %s", rhs_filename);
			fprintf(stderr, "[IO] Loading b from %s\n", rhs_filename);
			for (int i = 0; i < n; i++) {
				if (1 != fscanf(f_b, "%lg\n", &b[i]))
					errx(1, "parse error entry %d\n", i);
			}
			fclose(f_b);
		}
		else{
			for (int i = 0; i < n; i++)
				b[i] = PRF(i, seed);
		}
	}

	MPI_Bcast(x , n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(b , n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);


	/* solve Ax == b */
	cg_solve(A2, b, x, THRESHOLD, my_rank, nb_proc, scratch);

	if(my_rank == 0){
		/* Check result */
		if (safety_check) {
			y = scratch;
			sp_gemv(A, x, y);	// y = Ax
			for (int i = 0; i < n; i++)	// y = Ax - b
				y[i] -= b[i];
			fprintf(stderr, "[check] max error = %2.2e\n", norm(n, y));
		}

		/* Dump the solution vector */
		FILE *f_x = stdout;
		if (solution_filename != NULL) {
			f_x = fopen(solution_filename, "w");
			if (f_x == NULL)
				err(1, "cannot open solution file %s", solution_filename);
			fprintf(stderr, "[IO] writing solution to %s\n", solution_filename);
		}
		for (int i = 0; i < n; i++)
			fprintf(f_x, "%a\n", x[i]);
	}

	//MPI_Type_free(&mpi_csr_matrix);
	MPI_Finalize();
	return EXIT_SUCCESS;
}
