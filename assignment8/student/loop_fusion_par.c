#include <stdio.h>
void compute(unsigned long **a, unsigned long **b, unsigned long **c, unsigned long **d, int N, int num_threads) {

	// perform loop fusion to transform this loop and parallelize it with OpenMP
	fprintf(stdout,"b1N:%ld\n",b[1][N]);
	#pragma omp parallel for
	for (int i = 1; i < N; i++) {
		c[i][0] = (2* b[i][0] - 2* b[i][2]);
		for (int j = 1; j < N; j++) {
			a[i][j] = 2 * b[i][j];
			d[i][j] = a[i][j] * c[i][j];
			if(j < (N-2)){
				c[i][j] = a[i][j] - (2 * b[i][j + 2]);
			}
		}
		c[i][N-2] = a[i][N-2] - a[i][N];
	}


	// for (int i = 1; i < N; i++) {
	// 	c[i][0] = 2 * (b[i][0] - b[i][2]);
	// 	for (int j = 1; j < N-2; j++) {
	// 		c[i][j] = a[i][j] - 2 * b[i][j + 2]; // use B on right one
	// 	}
	// 	c[i][N-2] = a[i][N-2] - a[i][N]; // use B on right one
	// }
}


void compute_org(unsigned long **a, unsigned long **b, unsigned long **c, unsigned long **d, int N, int num_threads) {

	// perform loop fusion to transform this loop and parallelize it with OpenMP
	for (int i = 1; i < N; i++) {
		for (int j = 1; j < N; j++) {
			a[i][j] = 2 * b[i][j];
			d[i][j] = a[i][j] * c[i][j];
		}
	}

	for (int j = 1; j < N; j++) {
		for (int i = 1; i < N; i++) {
			c[i][j - 1] = a[i][j - 1] - a[i][j + 1];
		}
	}
}
