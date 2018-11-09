#ifndef DGEQP3_H
#define DGEQP3_H

#include "QR_utils.h"

void exec_dgeqp3(int m, int n, double *A, int lda, int *JPVT, double *TAU)
{
	int INFO;
    int LWORK = -1;
	
    // LWORK je -1 pa ovaj poziv dgeqp3 vraca samo optimalnu velicinu polja WORK u lwork2[0].
    double lwork2;    
    dgeqp3(&m, &n, A, &lda, JPVT, TAU, &lwork2, &LWORK, &INFO);

    int lwork3 = (int)lwork2;
    double *WORK = resize_WORK(lwork3);
    // izvrsavanje QRCP faktorizacije
    dgeqp3(&m, &n, A, &lda, JPVT, TAU, WORK, &lwork3, &INFO);
    if(INFO != 0)
	{
		fprintf(stderr,"Greska pri izvrsavanju QR faktorizacije! Kod greske je %d.\n",INFO);
		exit(1);
	}
}

#endif