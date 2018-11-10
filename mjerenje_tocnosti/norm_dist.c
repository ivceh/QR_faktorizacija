#include <stdio.h>
#include <stdlib.h>
#include "../QR_faktorizacije.h"
#include "mjerenje_tocnosti_utils.h"

/*
Argumenti komandne linije:
	argv[1]: m, visina matrice
	argv[2]: n, sirina matrice
	argv[3]: r, aproksimacijski rang
	argv[4] i argv[5]: b i p, parametri u randomiziranoj QRCP faktorizaciji
*/
int main(int argc, char *argv[])
{
	int i, j, m, n, r, b, p;
    double start, end;
	
	if(argc!=6)
	{
		fprintf(stderr, "Greska! Krivi broj argumenata komandne linije.\n");
		exit(-1);
	}
	
	m = atoi(argv[1]);
	n = atoi(argv[2]);
	r = atoi(argv[3]);
	b = atoi(argv[4]);
	p = atoi(argv[5]);
	double *A = calloc(m*n, sizeof(double)),
		   *B = calloc(m*n, sizeof(double));
	
	norm_dist_vect(A, m*n);

    int rowcolmin = (m<n) ? m : n;
    double TAU[rowcolmin];
	int JPVT[n];
	
	
	// dgeqrf
	copy_matrix(m, n, A, B);
	exec_dgeqrf(m, n, B, m, TAU);
	low_rank_approx_from_QR(m, n, r, B, m, TAU);
	printf("dgeqrf:  %g\n", F_distance(B, A, m, n));
	
	// dgeqp3
	copy_matrix(m, n, A, B);
	exec_dgeqp3(m, n, B, m, JPVT, TAU);
	low_rank_approx_from_QRCP(m, n, r, B, m, JPVT, TAU);
	printf("dgeqp3:  %g\n", F_distance(B, A, m, n));
	
	// rsrqrcp
	copy_matrix(m, n, A, B);
	rsrqrcp(m, n, B, JPVT, TAU, b, p);
	low_rank_approx_from_QRCP(m, n, r, B, m, JPVT, TAU);
	printf("rsrqrcp: %g\n", F_distance(B, A, m, n));
	
	// rqrcp
	copy_matrix(m, n, A, B);
	rqrcp(m, n, B, JPVT, TAU, b, p);
	low_rank_approx_from_QRCP(m, n, r, B, m, JPVT, TAU);
	printf("rqrcp:   %g\n", F_distance(B, A, m, n));
	
	// svd
	copy_matrix(m, n, A, B);
	low_rank_approx_from_SVD(m, n, r, B);
	printf("svd:     %g\n", F_distance(B, A, m, n));
	
	return 0;
}