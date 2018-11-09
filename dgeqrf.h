#ifndef DGEQRF_H
#define DGEQRF_H

#include "QR_utils.h"

void exec_dgeqrf(int m, int n, double *A, int LDA, double *TAU)
{
	int INFO;
    int LWORK = -1;
	
    // LWORK je -1 pa ovaj poziv dgeqrf vraca samo optimalnu velicinu polja WORK u lwork2[0].
    double lwork2;    
    dgeqrf(&m, &n, A, &LDA, TAU, &lwork2, &LWORK, &INFO);

    int lwork3 = (int)lwork2;
    double *WORK = resize_WORK(lwork3);
    // izvrsavanje QR faktorizacije
    dgeqrf(&m, &n, A, &LDA, TAU, WORK, &lwork3, &INFO);
    if(INFO != 0)
	{
		fprintf(stderr,"Greska pri izvrsavanju QR faktorizacije! Kod greske je %d.\n",INFO);
		exit(1);
	}
}

#endif