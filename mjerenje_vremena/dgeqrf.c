#include <stdio.h>
#include <stdlib.h>
#include "../QR_faktorizacije.h"

// argumenti komandne linije su dimenzije matrice A
int main(int argc, char *argv[])
{
    int i, j, m, n;
    double start, end;
	
	if(argc!=3)
	{
		fprintf(stderr, "Greska! Krivi broj argumenata komandne linije.\n");
		exit(-1);
	}
	
	m = atoi(argv[1]);
	n = atoi(argv[2]);
	double *A = calloc(m*n, sizeof(double));
	
	norm_dist_vect(A, m*n);

    int rowcolmin = (m<n) ? m : n;
    double TAU[rowcolmin];

	start = omp_get_wtime();
	
    exec_dgeqrf(m, n, A, m, TAU);
	
	end = omp_get_wtime();
	
	printf("Vrijeme: %lg s\n", end-start);
	
	free(A);
	resize_WORK(-1);
    return 0;
}