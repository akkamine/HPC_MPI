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
 *			$ ./cg --matrix bcsstk13.mtx                # loading matrix from file
 *			$ ./cg --matrix bcsstk13.mtx > /dev/null    # ignoring solution
 *			$ ./cg < bcsstk13.mtx > /dev/null           # loading matrix from stdin
 *      $  zcat matrix.mtx.gz | ./cg                # loading gziped matrix from
 *      $ ./cg --matrix bcsstk13.mtx --seed 42      # changing right-hand side
 *      $ ./cg --no-check < bcsstk13.mtx            # no safety check
 *
 * PRO-TIP :
 *			# downloading and uncompressing the matrix on the fly
 *			$ curl --silent https://hpc.fil.cool/matrix/bcsstk13.mtx.gz | zcat | ./cg
 */
#include <stdio.h>
#include <stdlib.h>
#include <err.h>
#include <math.h>
#include <getopt.h>
#include <sys/time.h>
#include <immintrin.h>
#include <mpi.h>

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
	#pragma omp parallel for
	for (int i = 0; i < n; i++) {
		d[i] = 0.0;
		for (int u = Ap[i]; u < Ap[i + 1]; u++){
			if (i == Aj[u]){
				d[i] += Ax[u];
			}
		}
	}
}

/* Matrix-vector product (with A in CSR format) : y = Ax */

void sp_gemv(const struct csr_matrix_t *A, const double *x, double *y, int n)
{
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
	double sum = 0.;
	#pragma omp parallel for simd reduction(+:sum)
	for (int i = 0; i < n; i++)
		sum += x[i] * y[i];
	return sum;
}

/* euclidean norm (a.k.a 2-norm) */

double norm(const int n, const double *x)
{
	return sqrt(dot(n, x, x));
}


/*********************** conjugate gradient algorithm *************************/

/* Solve Ax == b (the solution is written in x) */
void cg_solve(const struct csr_matrix_t *A, const double *b, double *x, const double epsilon, int num_of_proc, int my_rank)
{
	MPI_Status status;
	MPI_Request request;

	/* to distribute (parts of) A */
	int tag_line_pointers  = 1;
	int tag_column_indexes = 2;
	int tag_coefficients   = 3;

	/* to distribute dense vectors */
	int tag_p  = 4;
	int tag_d  = 5;
	int tag_r  = 6;
	int tag_x  = 7;


	int n  = 0; //initialize to silence warnings
	int nz = 0; //initialize to silence warnings
	/* share variables n and nz */
	if(my_rank == 0){
		fprintf(stderr, "%-50s", "[CG] Communicating variables n and nz... ");
		n = A->n;
		nz = A->nz;
	}
	MPI_Bcast(&n , 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&nz, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(my_rank == 0)
		fprintf(stderr, "OK\n");

	/* all processes keep a record of how many lines are allocated
	   to each process and the index of each process' starting line */

	int rem = n % num_of_proc; //first rem lines will be assigned 1 extra line
	int * table1; // will contain the number of lines assigned to each process
	int * table2; // will contain the index of the first line assigned to each process
	table1 = (int *)malloc(sizeof(int) * num_of_proc);
	table2 = (int *)malloc(sizeof(int) * num_of_proc);

	if(rem == 0){
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


	/*
	 * in this next section:
	 * distribute the corresponding number of lines of A to each process
	 * also adjust the values of the line pointers (Ap) because
	 * each of the processes will only have a subset of the lines
	 */

	int *Ap;
 	int *Aj;
 	double *Ax;

	int my_number_of_lines = table1[my_rank];
	int my_starting_index = table2[my_rank];
	int pointers_to_receive = my_number_of_lines +1; //+1 to compute number of coefficients of the last line
	int my_number_of_coeffs;
	int offset; //because of the csr_matrix_t line pointers
	int * table3 = NULL; //only used by process 0

	if(my_rank == 0){
		fprintf(stderr, "%-50s", "[CG] Communicating number of coefficients... ");
		Ap = A->Ap;
		Aj = A->Aj;
		Ax = A->Ax;
		table3 = (int *)malloc(sizeof(int) * num_of_proc);
		for(int i = 0; i < num_of_proc; i++){
			table3[i] = Ap[table2[i] + table1[i]] - Ap[table2[i]];
		}
	}
	/* each process receives its number of coefficients */
	MPI_Scatter(table3, 1, MPI_INT, &my_number_of_coeffs, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(my_rank == 0)
		fprintf(stderr, "OK\n");

	if(my_rank == 0)
		fprintf(stderr, "%-50s", "[CG] All processes allocating memory... ");
	else{
		/* allocate memory to store this process' assigned part of the martix */
		Ap = (int *)malloc(sizeof(int) * pointers_to_receive);
		Aj = (int *)malloc(sizeof(int) * my_number_of_coeffs);
		Ax = (double *)malloc(sizeof(double) * my_number_of_coeffs);
	}
	if(my_rank == 0)
		fprintf(stderr, "OK\n");



	if(my_rank == 0){
		fprintf(stderr, "%-50s", "[CG] Distributing matrix coefficients... ");
		/* process 0 will be charged with the first lines */
		for(int i = 1; i < num_of_proc; i++){
			int number_of_pointers_for_proc = table1[i] +1; //+1 in order to compute number of coefficients on the last line
			int number_of_coeffs_for_proc = table3[i];
			MPI_Ssend(&Ap[table2[i]],     number_of_pointers_for_proc, MPI_INT,    i, tag_line_pointers,  MPI_COMM_WORLD);
			MPI_Ssend(&Aj[Ap[table2[i]]], number_of_coeffs_for_proc,   MPI_INT,    i, tag_column_indexes, MPI_COMM_WORLD);
			MPI_Ssend(&Ax[Ap[table2[i]]], number_of_coeffs_for_proc,   MPI_DOUBLE, i, tag_coefficients,   MPI_COMM_WORLD);
		}

	}
	else{
		MPI_Recv(Ap, pointers_to_receive, MPI_INT,    0, tag_line_pointers,  MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(Aj, my_number_of_coeffs, MPI_INT,    0, tag_column_indexes, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(Ax, my_number_of_coeffs, MPI_DOUBLE, 0, tag_coefficients,   MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	/* indexes on the original matrix do not correspond to assigned lines,
		 pointer to the first assigned line therefore gives you the offset */
	offset = Ap[0];
	#pragma omp parallel for simd
	for(int i = 0; i < my_number_of_lines+1; i++)
		Ap[i] -= offset;

	if(my_rank == 0){
		free(table3);
		fprintf(stderr, "OK\n");
	}
	/* end of distribution of lines of A (+ columns and indexes) */


	/*
	 * in this next section:
	 * - all processes allocate memory for vectors p, q, x, r, d
	 * - process 0 extract the diagonal
	 * - process 0 sends all required vector sections to other processes
	 * - process 0 also computes rz and norm_r
	 */

	double *r;
 	double *p;
 	double *q;
 	double *d;
	int total_size;

	if(my_rank == 0)
		total_size = n * 2 + my_number_of_lines * 2;
	else
		total_size = n + my_number_of_lines * 4;
	/* scratch contains the necessary memory for all algorithm vectors */
	double * scratch = (double *)malloc(sizeof(double) * total_size);
	/*
	 * Explanation:
	 * ----------------------------------
	 * process 0 needs memory for:
	 * - full p vector (for the matrix-vector multiplication)
	 * - full d vector (because it extracts the diagonal)
	 * - personal part of vectors r and q
	 * personal part of these last vectors is of length: my_number_of_lines
	 * vectors x and b are already allocated
	 * ----------------------------------
	 * processes > 0 need memory for:
	 * - full p vector (for the matrix-vector multiplication)
	 * - personal part of vectors d, x, r and q
	 * personal part of these last vectors is of length: my_number_of_lines
	 */

	double rz;
	double norm_r;

	if(my_rank == 0){
		fprintf(stderr, "%-50s", "[CG] Distributing dense vectors... ");
		p = scratch; //full vector of length n
		d = scratch + n; //full vector of length n (for process 0 only)
		r = scratch + n + n;
		q = scratch + n + n + my_number_of_lines;
		double *z = scratch + n;

		extract_diagonal(A, d);
		/* norm_r is the scalar used as the
		   condition in the main while loop */
		rz = 0.;
		norm_r = 0.;
		#pragma omp parallel for simd reduction(+:rz,norm_r)
		for(int i = 0; i < n; i++){
			p[i] = b[i] / d[i]; //p <- r ./ d but at the beginning r = b, so this is equivalent
			double tmp_res = b[i] * b[i];
			norm_r += tmp_res; //norm(r) = norm(b) because r = b at the beginning
			rz += tmp_res / d[i]; //compute rz and norm_r now to avoid extra MPI communication
		}
		norm_r = sqrt(norm_r);

		/*
		 * r = b at the beginning of the algorithm
		 * only process 0 needs to do the following for loop
		 * processes > 0 will receive their subset of r directly (see below)
		 */
		#pragma omp parallel for simd
		for(int i = 0; i < my_number_of_lines; i++)
			r[i] = b[i];

		/*
		 * full vector p communicated to every process
		 * for vectors d and r, send only section that concerns them
		 * NOTE: at the beginning, r = b, so we send sections of b (equivalent to sections of r)
		 */
		for(int i = 1; i < num_of_proc; i++){
			int number_of_lines_for_proc = table1[i];
			MPI_Issend(p,             n,                        MPI_DOUBLE, i, tag_p, MPI_COMM_WORLD, &request);
			MPI_Issend(&d[table2[i]], number_of_lines_for_proc, MPI_DOUBLE, i, tag_d, MPI_COMM_WORLD, &request);
			MPI_Issend(&b[table2[i]], number_of_lines_for_proc, MPI_DOUBLE, i, tag_r, MPI_COMM_WORLD, &request);
		}
		for(int i = 1; i < num_of_proc; i++){
			MPI_Wait(&request, MPI_STATUS_IGNORE);
			MPI_Wait(&request, MPI_STATUS_IGNORE);
			MPI_Wait(&request, MPI_STATUS_IGNORE);
		}
		fprintf(stderr, "OK\n");
		/* end of section for process 0 */
	}
	else{
		p = scratch; //full vector of length n
		d = scratch + n;
		x = scratch + n + my_number_of_lines;
		r = scratch + n + my_number_of_lines * 2;
		q = scratch + n + my_number_of_lines * 3;
		MPI_Recv(p, n,                  MPI_DOUBLE, 0, tag_p, MPI_COMM_WORLD, &status);
		MPI_Recv(d, my_number_of_lines, MPI_DOUBLE, 0, tag_d, MPI_COMM_WORLD, &status);
		MPI_Recv(r, my_number_of_lines, MPI_DOUBLE, 0, tag_r, MPI_COMM_WORLD, &status);
		/* end of section for processes > 0 */
	}

	/* broadcast starting values of rz and norm_r */
	if(my_rank == 0)
		fprintf(stderr, "%-50s", "[CG] Broadcasting rz and norm(r)... ");
	MPI_Bcast(&norm_r, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&rz    , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if(my_rank == 0)
		fprintf(stderr, "OK\n");
	/* end of distribution of initial vectors and scalar values */



	/*
	 * This function follows closely the pseudo-code given in the (english)
	 * Wikipedia page "Conjugate gradient method". This is the version with
	 * preconditionning.
	 */

	/* We use x == 0 --- this avoids the first matrix-vector product. */
	#pragma omp parallel for simd
	for(int i = 0; i < my_number_of_lines; i++)
		x[i] = 0.;

	double * tmp; //used by process 0 for MPI_Gather
	if(my_rank == 0)
		tmp = (double *)malloc(sizeof(double) * num_of_proc * 2);
	else
		tmp = NULL;

	//scalars in one bloc to do Gather and Bcast more efficiently
	double scalar_bloc[3];
	//scalar_bloc[0] will be used to store rz
	//scalar_bloc[1] will be used to store norm_r
	//scalar_bloc[2] will be used to store beta

	struct csr_matrix_t A_sub; //to pass arguments to sp_gemv
	A_sub.n  = n;
	A_sub.nz = nz; //not really needed
	A_sub.Ap = Ap;
	A_sub.Aj = Aj;
	A_sub.Ax = Ax;

	if(my_rank == 0){
		fprintf(stderr, "[CG] Starting iterative solver\n");
		fprintf(stderr, "     ---> Working set : %.1fMbyte\n", 1e-6 * (12.0 * nz + 52.0 * n));
		fprintf(stderr, "     ---> Per iteration: %.2g FLOP in sp_gemv() and %.2g FLOP in the rest\n", 2. * nz, 12. * n);
	}
	double start = wtime();
	double last_display = start;
	int iter = 0;

	while(norm_r > epsilon){
		double old_rz = rz;
		/* each process calculates part of alpha */
		sp_gemv(&A_sub, p, q, my_number_of_lines); //ici
		double partial_pq = dot(my_number_of_lines, p+my_starting_index, q);

		/* process 0 computes and distributes alpha */
		MPI_Gather(&partial_pq, 1, MPI_DOUBLE, tmp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		double alpha;
		if(my_rank == 0){
			/* no parallel section as tmp is relatively short -> avoid unnecessary synchronizations */
			double pq = 0.;
			#pragma omp simd reduction(+:pq)
			for(int i = 0; i < num_of_proc; i++)
				pq += tmp[i];
			alpha = old_rz / pq;
		}
		MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		/* each process calculates part of rz and norm_r
		   and also updates personal subsection of x and r */
		rz = 0.;
		norm_r = 0.;
		#pragma omp simd reduction(+:rz,norm_r)
		for(int i = 0; i < my_number_of_lines; i++){
			x[i] += alpha * p[i + my_starting_index];
			r[i] -= alpha * q[i];
			double tmp_res = r[i] * r[i];
			norm_r += tmp_res;
			rz += tmp_res / d[i];
		}
		scalar_bloc[0] = rz;
		scalar_bloc[1] = norm_r;
		MPI_Gather(scalar_bloc, 2, MPI_DOUBLE, tmp, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		/* process 0 computes and distributes rz, norm_r and beta */
		double beta;
		if(my_rank == 0){
			/* no parallel section as tmp is relatively short -> avoid unnecessary synchronizations */
			rz = 0.;
			norm_r = 0.;
			#pragma omp simd reduction(+:rz,norm_r)
			for(int i = 0; i < num_of_proc*2; i+=2){
				rz += tmp[i];
				norm_r += tmp[i+1];
			}
			norm_r = sqrt(norm_r);
			beta = rz / old_rz;
			scalar_bloc[0] = rz;
			scalar_bloc[1] = norm_r;
			scalar_bloc[2] = beta;
		}

		/* transmit all 3 variables at once to reduce number of messages */
		MPI_Bcast(&scalar_bloc, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		rz     = scalar_bloc[0];
		norm_r = scalar_bloc[1];
		beta   = scalar_bloc[2];

		/* each process updates his section of p
		   and then broadcasts it */

		#pragma omp parallel for simd
		for(int i = 0; i < my_number_of_lines; i++)
			p[i+my_starting_index] = r[i] / d[i] + beta * p[i+my_starting_index];

		for(int i = 0; i < num_of_proc; i++)
			MPI_Bcast(&p[table2[i]], table1[i], MPI_DOUBLE, i, MPI_COMM_WORLD);

		iter++;
		double t = wtime();
		if(t - last_display > 0.5 && my_rank == 0){
			/* verbosity */
			double rate = iter / (t - start);	// iterations per s.
			double GFLOPs = 1e-9 * rate * (2 * nz + 12 * n);
			fprintf(stderr, "\r     ---> error : %2.2e, iter : %d (%.1f it/s, %.2f GFLOPs)", norm_r, iter, rate, GFLOPs);
			fflush(stdout);
			last_display = t;
		}
	} //end of main while loop

	/* assemble complete x vector */
	if(my_rank == 0){
		for(int i = 1; i < num_of_proc; i++)
			MPI_Recv(&x[table2[i]], table1[i], MPI_DOUBLE, i, tag_x, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		//for(int i = 1; i < num_of_proc; i++)
		//	MPI_Wait(&request, MPI_STATUS_IGNORE);
	}
	else{
		MPI_Ssend(x, my_number_of_lines, MPI_DOUBLE, 0, tag_x, MPI_COMM_WORLD);
	}

	free(table1);
	free(table2);
	free(scratch);

	if(my_rank == 0){
		free(tmp); //tmp == NULL for other processes, no need to free
		fprintf(stderr, "\n     ---> Finished in %.1fs and %d iterations\n", wtime() - start, iter);
	}
	else{
		free(Ap);
		free(Aj);
		free(Ax);
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

	int num_of_proc; // number of processes
	int my_rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_of_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	struct csr_matrix_t *A = NULL;
	double *x = NULL;
	double *b = NULL;
	double *y = NULL;
	int n = 0;

	/* Load the matrix */
	FILE *f_mat = stdin;
	if(my_rank == 0){
		if(matrix_filename){
			/* Load the matrix */
			f_mat = fopen(matrix_filename, "r");
			if(f_mat == NULL)
				err(1, "cannot matrix file %s", matrix_filename);
		}
		A = load_mm(f_mat);

		/* Allocate memory */
		n = A->n;
		double *mem = malloc(3 * n * sizeof(double));
		if (mem == NULL)
			err(1, "cannot allocate dense vectors");
		x = mem; /* solution vector */
		b = mem + n; /* right-hand side */
		y = mem + n + n; /* for safety check */
	}

	/* Prepare right-hand size */
	/* only process 0 will load/create vector b here */
	if(my_rank == 0){
		if(rhs_filename){	/* load from file */
			FILE *f_b = fopen(rhs_filename, "r");
			if(f_b == NULL)
				err(1, "cannot open %s", rhs_filename);
			fprintf(stderr, "[IO] Loading b from %s\n", rhs_filename);
			for(int i = 0; i < n; i++){
				if(1 != fscanf(f_b, "%lg\n", &b[i]))
					errx(1, "parse error entry %d\n", i);
			}
			fclose(f_b);
		}
		else{
			#pragma omp parallel for
			for (int i = 0; i < n; i++)
				b[i] = PRF(i, seed);
		}
	}

	/* solve Ax == b */
	cg_solve(A, b, x, THRESHOLD, num_of_proc, my_rank);

	/* Check result */
	if(my_rank == 0){
		/* done by process 0 for simplicity */
		if(safety_check){
			sp_gemv(A, x, y, n);	// y = Ax //ici
			#pragma omp parallel for simd
			for (int i = 0; i < n; i++)	// y = Ax - b
				y[i] -= b[i];
			fprintf(stderr, "[check] max error = %2.2e\n", norm(n, y));
		}
		/* Dump the solution vector */
		FILE *f_x = stdout;
		if(solution_filename != NULL){
			f_x = fopen(solution_filename, "w");
			if (f_x == NULL)
				err(1, "cannot open solution file %s", solution_filename);
			fprintf(stderr, "[IO] writing solution to %s\n", solution_filename);
		}
		for(int i = 0; i < n; i++)
			fprintf(f_x, "%a\n", x[i]);
	}

	MPI_Finalize();
	return EXIT_SUCCESS;
}
