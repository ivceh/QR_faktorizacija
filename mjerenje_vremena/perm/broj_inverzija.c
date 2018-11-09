#include <stdio.h>
#include <stdlib.h>

// racuna broj inverzija P i sortira P merge sortom
int broj_inverzija(int P[], int n)
{
	int kopija[n], brojac = 0, i, l, d;

	if(n<=1)
		return 0;
	else
	{
		for (i=0; i<n; ++i)
			kopija[i] = P[i];

		// pozivanje rekurzije na podnizovima
		brojac += broj_inverzija(kopija, n/2);
		brojac += broj_inverzija(kopija+n/2, (n+1)/2);

		// spajanje
		l = 0;
		d = n/2;
		i = 0;
		while (l<n/2 && d<n)
		{
			if (kopija[l] < kopija[d])
				P[i++] = kopija[l++];
			else
			{
				brojac += n/2 - l;
				P[i++] = kopija[d++];
			}
		}
		while (l<n/2)
			P[i++] = kopija[l++];
		while (d<n)
			P[i++] = kopija[d++];

		return brojac;
	}
}

// P se mijenja pozivom funkcije
void inverz(int P[], int n)
{
	int Inv[n];

	for (int i=0; i<n; ++i)
		Inv[P[i]-1] = i;

	for (int i=0; i<n; ++i)
		P[i] = Inv[i];
}

// Q = P Q (Q se mijenja pozivom funkcije)
void kompozicija(int P[], int Q[], int n)
{
	for (int i=0; i<n; ++i)
		Q[i] = P[Q[i]-1];
}

int main(int argc, char *argv[])
{
	int *P, *Q, n, i=0, r, maxr;

	if (argc != 2 && argc != 3)
	{
		fprintf(stderr, "Greska! Krivi broj argumenata komandne linije.");
		exit(-1);
	}

	FILE *fp = fopen(argv[1], "r");
	while(fscanf(fp, "%d", &n) == 1)
		++i;
	n = i;
	fseek(fp, 0, SEEK_SET);
	P = calloc(n, sizeof(int));
	for (i=0; i<n; ++i)
		fscanf(fp, "%d", P+i);
	fclose(fp);

	if (argc == 3)
	{
		fopen(argv[2], "r");
		Q = calloc(n, sizeof(int));
		for (i=0; i<n; ++i)
			fscanf(fp, "%d", Q+i);
		fclose(fp);

		inverz(P, n);
		kompozicija(P, Q, n);
		r = broj_inverzija(Q, n);
		maxr = n*(n-1)/2;
		printf("Inverzijska razlicitost je %d / %d = %g%%\n", r, maxr, (double)r/maxr*100);
	}
	else
	{
	    r = broj_inverzija(P, n);
		maxr = n*(n-1)/2;
		printf("Broj inverzija je %d / %d = %g%%\n", r, maxr, (double)r/maxr*100);
	}

	return 0;
}
