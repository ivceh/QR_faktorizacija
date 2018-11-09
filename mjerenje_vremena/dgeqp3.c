#include <stdio.h>
#include <stdlib.h>
#include "../QR_faktorizacije.h"

// argumenti komandne linije su dimenzije matrice A i opcionalno ime datoteke u koju upisujemo vektor permutacije JPVT
int main(int argc, char *argv[])
{
    int i, j, m, n;
    double start, end;
	
	if(argc!=3 && argc!=4)
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
	int JPVT[n];

	start = omp_get_wtime();
	
    exec_dgeqp3(m, n, A, m, JPVT, TAU);
	
	end = omp_get_wtime();
	
	if(argc == 4)
	{
		FILE *fp = fopen(argv[3], "w");
		if (fp == NULL)
		{
			fprintf(stderr, "Greska! Ne mogu otvoriti datoteku \"%s\" za pisanje.\n", argv[3]);
			exit(-1);
		}
		for (i=0; i<n; ++i)
			fprintf(fp, "%d ", JPVT[i]);
		fclose(fp);
	}
	
	printf("Vrijeme: %lg s\n", end-start);
	
	free(A);
	resize_WORK(-1);
    return 0;
}