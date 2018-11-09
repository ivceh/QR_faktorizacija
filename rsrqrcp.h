#ifndef RSRQRCP_H
#define RSRQRCP_H

#include "QR_utils.h"
#include "dgeqrf.h"
#include "dgeqp3.h"

double rsrqrcp(int m, int n, double *A, int *JPVT, double *TAU, int b, int p)
{
	int i, j, w, *JPVT2 = calloc(n, sizeof(int)), rowcolmin;
	double *Omega = resize_WORK((b+p)*m), *B = calloc((b+p)*n, sizeof(double)), *TAU2;
	
	for (i=0; i<n; ++i)
		JPVT[i] = i+1;
	
	for(i=0; i<n; i+=b)
	{
		w = min(n-i, b);
		
		if(w > 1 && b+p<m-i)
		{
			// generiraj slucajnu matricu dimenzija w+l, m-i
			Omega = resize_WORK((w+p)*(m-i));
			norm_dist_vect(Omega, (w+p)*(m-i));
			
			// B = Omega * A
			matr_mult(w+p, n-i, m-i, Omega, w+p, A+(m+1)*i, m, B, w+p);
			
			// QRCP faktorizacija matrice B
			for (j=0; j<n-i; ++j)
				JPVT2[j] = 0;
			rowcolmin = min(w+p, n-i);
			TAU2 = resize_TAU2(rowcolmin);
			exec_dgeqp3(w+p, n-i, B, w+p, JPVT2, TAU2);
			
			// permutiraj stupce matrice A prema permutaciji JPVT2
			permute_columns(m, n-i, A+i*m, JPVT2);
			
			// azuriraj permutaciju JPVT koristeci JPVT2
			update_perm(JPVT+i, n-i, JPVT2);
			
			// izvrsi QR faktorizaciju sljedecih w stupaca u A
			exec_dgeqrf(m-i, w, A+(m+1)*i, m, TAU+i);
			
			// azuriraj A desno od tih stupaca (po potrebi)
			if(i+b<n)
			{
				exec_dormqr(m-i, n-i-b, b, A+(m+1)*i, m, TAU+i, A+(i+b)*m+i, m);
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