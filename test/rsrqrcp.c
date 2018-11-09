#include <stdio.h>
#include <stdlib.h>
#include "../QR_faktorizacije.h"

/* 
   argumenti komandne linije:
   argv[1] = ime tekstualne datoteke u koju je spremljena matrica A, prvo pišu broj redaka i broj stupaca, a zatim elementi
   argv[2] = ime tekstualne datoteke u koju zelimo spremiti faktoriziranu matricu A (u gornjem trokutu R, a ispod njega vektori Householderovih reflektora)
   argv[3] = ime tekstualne datoteke u koju zelimo spremiti vektor permutacije JPVT
   argv[4] = ime tekstualne datoteke u koju zelimo spremiti TAU (vektor koeficijenata Householderovih reflektora)
   argv[5] = ime tekstualne datoteke u koju zelimo spremiti matricu A rekonstruiranu iz njene QR faktorizacije
   argv[6] = parametar b
   argv[7] = parametar p
*/

int main(int argc, char *argv[])
{
	FILE *fp;
    int i, j, m, n, b, p;
	
	if(argc!=8)
	{
		fprintf(stderr, "Greska! Krivi broj argumenata komandne linije.\n");
		exit(-1);
	}
	
	b = atoi(argv[6]);
	p = atoi(argv[7]);
	
	fp = fopen(argv[1], "r");
	if (fp == NULL)
	{
		fprintf(stderr, "Greska! Ne mogu otvoriti datoteku \"%s\" za citanje.\n", argv[1]);
		exit(-1);
	}
	// ucitavanje dimenzija matrice A
	if (fscanf(fp, "%d %d", &m, &n) != 2)
	{
		fprintf(stderr, "Greska pri ucitavanju dimenzija matrrice!\n");
		exit(-1);
	}
	double A[m*n];
    // ucitavanje matrice A
	for (i=0; i<m; ++i)
		for (j=0; j<n; ++j)
			if (fscanf(fp, "%lg", A+(j*m+i)) != 1)
			{
				fprintf(stderr, "Greska pri ucitavanju matrice!\n");
				exit(-1);
			}
	fclose(fp);

	
    int rowcolmin = min(m, n);
    double TAU[rowcolmin];
	int JPVT[n];
	
	// izvodenje algoritma
    rsrqrcp(m, n, A, JPVT, TAU, b, p);
	
	
	fp = fopen(argv[2], "w");
	if (fp == NULL)
	{
		fprintf(stderr, "Greska! Ne mogu otvoriti datoteku \"%s\" za pisanje.\n", argv[2]);
		exit(-1);
	}
	// ispis matrice A: u gornjem trokutu R a ispod njega vektori Householderovih reflektora
	for(i=0; i<m; ++i)
	{
		for(j=0; j<n; ++j)
			fprintf(fp, "%7g ", A[j*m+i]);
		fprintf(fp, "\n");
	}
	fclose(fp);
	
	fp = fopen(argv[3], "w");
	if (fp == NULL)
	{
		fprintf(stderr, "Greska! Ne mogu otvoriti datoteku \"%s\" za pisanje.\n", argv[3]);
		exit(-1);
	}
	// ispis vektora JPVT
	for(i=0; i<n; ++i)
		fprintf(fp, "%d ", JPVT[i]);
	fclose(fp);
	
	fp = fopen(argv[4], "w");
	if (fp == NULL)
	{
		fprintf(stderr, "Greska! Ne mogu otvoriti datoteku \"%s\" za pisanje.\n", argv[4]);
		exit(-1);
	}
	// ispis vektora TAU
	for(i=0; i<rowcolmin; ++i)
		fprintf(fp, "%g ", TAU[i]);
	fclose(fp);
	
	
	// rekonstrukcija matrice A
	construct_A_from_QRCP(m, n, A, m, JPVT, TAU);
	
	fp = fopen(argv[5], "w");
	if (fp == NULL)
	{
		fprintf(stderr, "Greska! Ne mogu otvoriti datoteku \"%s\" za pisanje.\n", argv[5]);
		exit(-1);
	}
	for (i=0; i<m; ++i)
	{
		for(j=0; j<n; ++j)
			fprintf(fp, "%7g ", A[j*m+i]);
		fprintf(fp, "\n");
	}
	fclose(fp);
	
    resize_WORK(-1);
	resize_TAU2(-1);
    return 0;
}