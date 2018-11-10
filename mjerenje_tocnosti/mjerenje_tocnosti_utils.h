#ifndef MJERENJE_TOCNOSTI_UTILS_H
#define MJERENJE_TOCNOSTI_UTILS_H

#include <mkl.h>
#include <omp.h>

double F_norm(double *A, int m, int n)
{
	char normF = 'F';
	
	return dlange(&normF, &m, &n, A, &m, NULL);
}

// A se mijenja pozivom ove funkcije!
double F_distance(double *A, double *B, int m, int n)
{
	#pragma omp parallel for
	for (int i=0; i<m*n; ++i)
		A[i] -= B[i];
	
	return F_norm(A, m, n);
}

void low_rank_approx_from_QR(int m, int n, int r, double *A, int lda, double *TAU)
{
	int INFO, LWORK = -1, k;
	char side='L', notransp='N';
	double *A_upper = calloc(m*n, sizeof(double));
	
	#pragma omp parallel for
	for (int k=0; k<m*n; ++k)
	{
		int i = k%m, j = k/m;
		if(i<=j && i<r)
			A_upper[k] = A[j*lda+i];
		else
			A_upper[k] = 0;
	}
	
	k = min(min(m, n), r);
	
    // LWORK je -1 pa ovaj poziv dormqr vraca samo optimalnu velicinu polja WORK u lwork2[0].
    double lwork2;    
    dormqr(&side, &notransp, &m, &n, &k, A, &lda, TAU, A_upper, &m, &lwork2, &LWORK, &INFO);

    int lwork3 = (int)lwork2;
    double *WORK = resize_WORK(lwork3);
    // mnozenje matricom Q^T slijeva
    dormqr(&side, &notransp, &m, &n, &k, A, &lda, TAU, A_upper, &m, WORK, &lwork3, &INFO);
    if(INFO != 0)
	{
		fprintf(stderr,"Greska pri mnozenju matricom Q^T! Kod greske je %d.\n",INFO);
		exit(1);
	}
	
	#pragma omp parallel for
	for (int k=0; k<m*n; ++k)
	{
		int i = k%m, j = k/m;
		A[j*lda+i] = A_upper[k];
	}
}

void low_rank_approx_from_QRCP(int m, int n, int r, double *A, int lda, int *JPVT, double *TAU)
{
	low_rank_approx_from_QR(m, n, r, A, lda, TAU);
	inv_permute_columns(m, n, A, JPVT);
}

void copy_matrix(int m, int n, double *A, double *B)
{
	#pragma omp parallel for
	for (int i=0; i<m*n; ++i)
		B[i] = A[i];
}

void predetermined_singular_values_matrix(int m, int n, double *A)
{
	#pragma omp parallel for
	for (int k=0; k<m*n; ++k)
	{
		int i = k%m, j = k/m;
		if(i == j)
		{
			// ovdje odredujemo singularne vrijednosti
			A[k] = (double)i/m;
		}
		else
			A[k] = 0;
	}
	
	static int seed[4] = {1,1,1,1};
	char sidel = 'L', sider = 'R', noinit = 'N';
	double *WORK = resize_WORK(3 * max(m,n));
	int INFO;
	
	dlaror(&sidel, &noinit, &m, &n, A, &m, seed, WORK, &INFO);
	if (INFO != 0)
	{
		fprintf(stderr, "Greska pri mnozenju slucajnom ortogonalnom matricom slijeva! Kod greske je %d\n", INFO);
		exit(-1);
	}
	
	dlaror(&sider, &noinit, &m, &n, A, &m, seed, WORK, &INFO);
	if (INFO != 0)
	{
		fprintf(stderr, "Greska pri mnozenju slucajnom ortogonalnom matricom slijeva! Kod greske je %d\n", INFO);
		exit(-1);
	}
}

void low_rank_approx_from_SVD(int m, int n, int r, double *A)
{
	char job = 'S', notransp = 'N', transp = 'T';
	int rowcolmin = min(m,n);
	double *u = calloc(m*rowcolmin, sizeof(double)),
		   *vt = calloc(n*rowcolmin, sizeof(double)),
		   S[rowcolmin];
	
	int INFO;
    int LWORK = -1;
	
	r = min(r, min(m, n));
	
    // LWORK je -1 pa ovaj poziv dgesvd vraca samo optimalnu velicinu polja WORK u lwork2[0].
    double lwork2;   
	dgesvd(&job, &job, &m, &n, A, &m, S, u, &m, vt, &rowcolmin, &lwork2, &LWORK, &INFO);

    int lwork3 = (int)lwork2;
    double *WORK = resize_WORK(lwork3);
    // izvrsavanje SVD faktorizacije
	dgesvd(&job, &job, &m, &n, A, &m, S, u, &m, vt, &rowcolmin, WORK, &lwork3, &INFO);
    if(INFO != 0)
	{
		fprintf(stderr,"Greska pri izvrsavanju QR faktorizacije! Kod greske je %d.\n",INFO);
		exit(1);
	}
	
	#pragma omp parallel for
	for (int i=0; i<m*n; ++i)
		A[i] = 0;
	
	#pragma omp parallel for
	for (int i=0; i<m*r; ++i)
		A[i] = u[i];
	free(u);
	
	#pragma omp parallel for
	for (int i=0; i<m*r; ++i)
		A[i] *= S[i/m];
	
	double *C = calloc(m*n, sizeof(double)), alpha = 1, beta = 0;
	
	dgemm(&notransp, &notransp, &m, &n, &r, &alpha, A, &m, vt, &rowcolmin, &beta, C, &m);
	
	#pragma omp parallel for
	for (int i=0; i<m*n; ++i)
		A[i] = C[i];
	
	free(C);
	free(vt);
}

#endif