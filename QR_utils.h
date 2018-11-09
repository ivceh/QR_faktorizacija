#ifndef QR_UTILS_H
#define QR_UTILS_H

#include <stdlib.h>
#include <mkl.h>
#include <omp.h>

int max(int a, int b)
{
	return (a>b) ? a : b;
}

int min(int a, int b)
{
	return (a<b) ? a : b;
}

double *resize_WORK(int n)
{
	static double *WORK = NULL;
	static int l = 0;
	
	if(n==-1)
	{
		free(WORK);
		return NULL;
	}
	else
	{
		if(n > l)
		{
			free(WORK);
			WORK = calloc(n, sizeof(double));
			l = n;
		}
		return WORK;
	}
}

double *resize_TAU2(int n)
{
	static double *TAU2 = NULL;
	static int l = 0;
	
	if(n==-1)
	{
		free(TAU2);
		return NULL;
	}
	else
	{
		if(n > l)
		{
			free(TAU2);
			TAU2 = calloc(n, sizeof(double));
			l = n;
		}
		return TAU2;
	}
}

void exec_dormqr(int m, int n, int k, double *A, int lda, double *TAU, double *C, int ldc)
{
	int INFO, LWORK = -1;
	char side='L', transp='T';
	
    // LWORK je -1 pa ovaj poziv dormqr vraca samo optimalnu velicinu polja WORK u lwork2[0].
    double lwork2;    
    dormqr(&side, &transp, &m, &n, &k, A, &lda, TAU, C, &ldc, &lwork2, &LWORK, &INFO);

    int lwork3 = (int)lwork2;
    double *WORK = resize_WORK(lwork3);
    // mnozenje matricom Q^T slijeva
    dormqr(&side, &transp, &m, &n, &k, A, &lda, TAU, C, &ldc, WORK, &lwork3, &INFO);
    if(INFO != 0)
	{
		fprintf(stderr,"Greska pri mnozenju matricom Q^T! Kod greske je %d.\n",INFO);
		exit(1);
	}
}

void norm_dist_vect(double *x, int n)
{
	int dist = 3, max_threads = omp_get_max_threads();
	static int *iseed, first_time = 1;
	
	if (first_time)
	{
		iseed = calloc(max_threads * 4, sizeof(int));
		for (int i=0; i<max_threads*4; ++i)
		{
			if (i%4 == 3)
				iseed[i] = i/2;
			else
				iseed[i] = 0;
		}
	}
	
	#pragma omp parallel
	{
		int i = omp_get_thread_num(),
			threads = omp_get_num_threads(),
			start, end, t, l;
		t = n / threads;
		start = i*t + min(i, n % threads);
		end = (i+1)*t + min(i+1, n % threads);
		l = end - start;
		dlarnv(&dist, iseed + i*4, &l, x+start); 
	}
	
	first_time = 0;
}

void matr_mult(int m, int n, int k, double *A, int lda, double *B, int ldb, double *C, int ldc)
{
	char notransp = 'N';
	double alpha = 1, beta = 0;
	dgemm(&notransp, &notransp, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
}

void permute_columns(int m, int n, double *A, int *JPVT)
{
	double *WORK = resize_WORK(m*n);
	
	#pragma omp parallel for
	for(int i=0; i<n; ++i)
	{
		int jpvti = JPVT[i]-1;
		for(int j=0; j<m; ++j)
			WORK[i*m+j] = A[jpvti*m+j];
	}
	
	#pragma omp parallel for
	for (int i=0; i<m*n; ++i)
		A[i] = WORK[i];
}

void inv_permute_columns(int m, int n, double *A, int *JPVT)
{
	double *WORK = resize_WORK(m*n);
	
	#pragma omp parallel for
	for(int i=0; i<n; ++i)
	{
		int jpvti = JPVT[i]-1;
		for(int j=0; j<m; ++j)
			WORK[jpvti*m+j] = A[i*m+j];
	}
	
	#pragma omp parallel for
	for (int i=0; i<m*n; ++i)
		A[i] = WORK[i];
}

void update_perm(int *JPVT, int l, int *JPVT2)
{
	int j;
	
	for (j=0; j<l; ++j)
		JPVT2[j] = JPVT[JPVT2[j]-1];
	for (j=0; j<l; ++j)
		JPVT[j] = JPVT2[j];
}

void construct_A_from_QR(int m, int n, double *A, int lda, double *TAU)
{
	int INFO, LWORK = -1, k=min(m,n);
	char side='L', notransp='N';
	double *A_upper = calloc(m*n, sizeof(double));
	
	#pragma omp parallel for
	for (int k=0; k<m*n; ++k)
	{
		int i = k%m, j = k/m;
		if(i<=j)
			A_upper[k] = A[j*lda+i];
		else
			A_upper[k] = 0;
	}
	
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

void construct_A_from_QRCP(int m, int n, double *A, int lda, int *JPVT, double *TAU)
{
	construct_A_from_QR(m, n, A, lda, TAU);
	inv_permute_columns(m, n, A, JPVT);
}

#endif