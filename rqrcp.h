#ifndef RQRCP_H
#define RQRCP_H

#include "QR_utils.h"
#include "dgeqrf.h"
#include "dgeqp3.h"

void update_B(int h, int b, int n, double *B, double *A, int lda)
{
	double *WORK = resize_WORK(b*n), one=1;
	int i, nminusb = n-b, *IPIV = calloc(b, sizeof(int)), INFO;
	char notransp='N';
	
	#pragma omp parallel for
	for (i=0; i<b*n; ++i)
		WORK[i] = A[(i/b)*lda+(i%b)];
	
	dgesv(&b, &nminusb, WORK, &b, IPIV, WORK+b*b, &b, &INFO);
	if(INFO != 0)
	{
		fprintf(stderr, "Greska pri rjesavanju sustava! Kod greske je %d.", INFO);
		exit(-1);
	}
	
	dgemm(&notransp, &notransp, &b, &nminusb, &b, &one, B, &h, WORK+b*b, &b, &one, B+b*h, &h);
	
	free(IPIV);
}

double rqrcp(int m, int n, double *A, int *JPVT, double *TAU, int b, int p)
{
	int i, j, w, *JPVT2 = calloc(n, sizeof(int)), rowcolmin;
	double *Omega = resize_WORK((b+p)*m), *B = calloc((b+p)*n, sizeof(double)), *TAU2;
	
	for (i=0; i<n; ++i)
		JPVT[i] = i+1;
	
	if (n < b)
		b = n;
	
	// generiraj slucajnu matricu dimenzija b+p, m
	norm_dist_vect(Omega, (b+p)*m);
	
	// B = Omega * A
	matr_mult(b+p, n, m, Omega, b+p, A, m, B, b+p);
	
	for(i=0; i<n; i+=b)
	{
		w = min(n-i, b);
		
		if(w > 1 && b+p<m-i)
		{
			// QRCP faktorizacija matrice B
			for (j=0; j<n-i; ++j)
				JPVT2[j] = 0;
			rowcolmin = min(w+p, n-i);
			TAU2 = resize_TAU2(rowcolmin);
			exec_dgeqp3(w+p, n-i, B+(b+p)*i, w+p, JPVT2, TAU2);
			
			// permutiraj stupce matrice A prema permutaciji JPVT2
			permute_columns(m, n-i, A+i*m, JPVT2);
			
			// azuriraj permutaciju JPVT koristeci JPVT2
			update_perm(JPVT+i, n-i, JPVT2);
			
			// izvrsi QR faktorizaciju sljedecih w stupaca u A
			exec_dgeqrf(m-i, w, A+(m+1)*i, m, TAU+i);
			
			if(i+b<n)
			{
				// azuriraj A desno od tih stupaca (po potrebi)
				exec_dormqr(m-i, n-i-b, b, A+(m+1)*i, m, TAU+i, A+(i+b)*m+i, m);
				
				// pripremi B za sljedeci korak
				update_B(b+p, b, n-i, B+(b+p)*i, A, m);
			}
		}
		else //randomizacija vise nema smisla jer se njome nista ne optimizira, provodimo obicnu QRCP faktorizaciju do kraja
		{
			// izvrsi QRCP faktorizaciju ostatka matrice A
			exec_dgeqp3(m-i, n-i, A+(m+1)*i, m, JPVT2, TAU+i);
			
			// azuriraj permutaciju JPVT koristeci JPVT2
			update_perm(JPVT+i, n-i, JPVT2);
			
			break;
		}
	}
	
	free(JPVT2);
	free(B);
}

#endif