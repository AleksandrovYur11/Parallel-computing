#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include "omp.h"

void merge(double M1[], double M2[], int sizeArrayM2) {
	for (int k = 0; k < sizeArrayM2; k++) {
		M2[k] = fabsf(M1[k] - M2[k]);
	}
}


int isSorted(double M2[], int sizeArrayM2)
{
	while (--sizeArrayM2 > 1)
		if (M2[sizeArrayM2] < M2[sizeArrayM2 - 1])
			return 0;
	return 1;
}

void swap(double* x, double* y)
{
	double tmp = *x;
	*x = *y;
	*y = tmp;
}

void shuffle(double M2[], int sizeArrayM2)
{
	for (int i = 0; i < sizeArrayM2; i++)
		swap(&M2[i], &M2[rand() % sizeArrayM2]);
}

void stupid_sort(double M2[], int sizeArrayM2)
{
	while (!isSorted(M2, sizeArrayM2))
		shuffle(M2, sizeArrayM2);
}

int main(int argc, char* argv[]) {

	const int A = 440;
	int numberFlow;
	int N;
	struct timeval T1, T2;
	long delta_ms;

	if (argc < 3) {
		printf("Need to enter arguments (sizeArray) and number flow ");
		return -1;
	}

	N = atoi(argv[1]);
	numberFlow = atoi(argv[2]);

#ifdef _OPENMP
	omp_set_num_threads(numberFlow);
#endif

	gettimeofday(&T1, NULL);

	unsigned int first_random = 0;
	unsigned int second_random = 0;

	for (int i = 0; i < 100; i++) {

		srand(i);


		double* M1 = (double*)malloc(N * sizeof(double));
		double* M2 = (double*)malloc(N / 2 * sizeof(double));
		double* M2_Copy = (double*)malloc(N / 2 * sizeof(double));


		for (int k = 0; k < N; k++) {
			double value = 1 + rand_r(&first_random) % A;
			M1[k] = value;
		};

		for (int k = 0; k < N / 2; k++) {
			double value = 1 + rand_r(&second_random) % (A * 10);
			M2[k] = value;
			M2_Copy[k] = value;
		}



#pragma omp parallel default(none)  shared(N, M1, M2, M2_Copy) 
		{

#pragma omp for  schedule(runtime)
			for (int k = 0; k < N; k++) {
				M1[k] = cosh(sqrt(M1[k])) / sinh(sqrt(M1[k]));
			}


#pragma omp for  schedule(runtime)
			for (int k = 0; k < N / 2; k++) {
				if (k == 0) {
					M2[k] = fabsf(cos(M2[k]));
					continue;
				}
				M2[k] = fabsf(cos(M2[k] + M2_Copy[k - 1]));
			}

#pragma omp for  schedule(runtime)
			for (int k = 0; k < N / 2; k++) {
				M2[k] = fabsf(M1[k] - M2[k]);
			}
		}

			stupid_sort(M2, N / 2);


			double result = 0;
			double min = 0;
			for (int i = 0; i < N / 2; i++) {
				if (M2[i] != 0) {
					min = M2[i];
					break;
				}
			}

#pragma omp parallel for default(none) shared(N, M2, min)  reduction(+:result) schedule(guided, 8)
			for (int k = 0; k < N / 2; k++) {
				if (((long)(M2[k] / min) % 2) == 0) {
					result += sin(M2[k]);
				}
			}

			printf("Result: %f\n", result);

			free(M1);
			free(M2);
			free(M2_Copy);
	}

		gettimeofday(&T2, NULL);
		delta_ms = 1000 * (T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec) / 1000;
		printf("\nN=%d. Milliseconds passed: %ld\n", N, delta_ms);

		return 0;
	
}


