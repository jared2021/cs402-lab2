#include "cblas.h"
#include <cstdio>

void compare_matrix(int n, double* A1, double* A2){
	int i, j;
	double d1, d2;
	for (i=0; i<n;i++){
		for(j=0;j<n;j++){
			d1= *(A1 + i*n + j);
			d2= *(A2 + i*n + j);
			if(fabs(d2-d1)/(d2) > 1e-6){
				printf("ERROR: %f(my)<>%f(dgemm)\n", d1, d2);
				exit(1);
			}
		}
	}
	printf("Correct result! :-) \n");
}

void init_matrix(int n, double* A){
	int i,j;
	for (i=0;i<n;i++){
		for(j=0;j<n;j++){
			*(A + i*n +j)= rand() / (RAND_MAX * 1.0);
		}
	}
}

void unoptimized_dgemm(int n, double* A, double* B, double* C)
{
	for (int i=0; i<n;i++)
	{
		for (int j=0; j<n;j++)
		{
			double cij= C[i+j*n]; /* cij= C[i][j] */
			for (int k=0;k<n;k++)
			{
				cij += A[i+k*n] * B[k+j*n]; /* cij += A[i][k] * B[k][j] */
			}
			C[i+j*n] = cij; /*C[i][j]= cij */
		}
	}
}
