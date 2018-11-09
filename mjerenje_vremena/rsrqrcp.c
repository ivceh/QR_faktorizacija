#include <stdio.h>
#include <stdlib.h>
#include "../QR_faktorizacije.h"

// argumenti komandne linije su dimenzije matrice A, b, p i opcionalno ime datoteke u koju upisujemo vektor permutacije JPVT
int main(int argc, char *argv[])
{
    int i, j, m, n, b, p;
    double start, end;
	
	if(argc!=5 && argc!=6)
	{
		fprintf(stderr, "Greska! Krivi broj argumenata komandne linije.\n");
		exit(-1);
	}
	m = atoi(argv[1]);
	n = atoi(argv[2]);
	b = atoi(argv[3]);
	p = atoi(argv[4]);
	double *A = calloc(m*n, sizeof(double));
	
	norm_dist_vect(A, m*n);

    int rowcolmin = (m<n) ? m : n;
    double TAU[rowcolmin];
	int JPVT[n];

	start = omp_get_wtime();
	
    rsrqrcp(m, n, A, JPVT, TAU, b, p);
	
	end = omp_get_wtime();
	
	if(argc == 6)
	{
		FILE *fp = fopen(argv[5], "w");
		if (fp == NULL)
		{
			fprintf(stderr, "Greska! Ne mogu otvoriti datoteku \"%s\" za pisanje.\n", argv[5]);
			exit(-1);
		}
		for (i=0; i<n; ++i)
			fprintf(fp, "%d ", JPVT[i]);
		fclose(fp);
	}
	
	printf("Vrijeme: %lg s\n", end-start);
	
	free(A);
	resize_WORK(-1);
	resize_TAU2(-1);
    return 0;
}