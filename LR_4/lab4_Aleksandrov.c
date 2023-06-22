#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include <float.h>
#ifdef _OPENMP
#include <omp.h>
#endif

void merge(double*  restrict M1, double* restrict M2, int sizeArrayM2) {
	for (int k = 0; k < sizeArrayM2; k++) {
		M2[k] = fabsf(M1[k] - M2[k]);
	}
}

int isSorted(double* restrict M2, int sizeArrayM2)
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

unsigned int f(unsigned int arg) {
	return ((int)&arg >> 3) * sqrt(arg) + 23;
}

void merge_sorted(double* restrict arr1, int n1, double* restrict arr2, int n2, double* restrict dst) {
	int i = 0, i1 = 0, i2 = 0;
	while (i < n1 + n2) {
		dst[i++] = arr1[i1] > arr2[i2] && i2 < n2 ? arr2[i2++] : arr1[i1++];
	}

}

void sort(double* restrict arr, int size, double* restrict dst) {
	
	
	int half = (int)size / 2;
#pragma omp parallel sections 
	{
#pragma omp section
		for (int i = 0; i < half; i++) {
			for (int j = 0; j < half; j++)
				if (arr[j] > arr[j + 1])
					swap(&arr[j], &arr[j + 1]);
		}
		
#pragma omp section
		for (int i = half; i < size - 1; i++) 
			for (int j = half; j < size - 1; j++)
				if (arr[j] > arr[j + 1])
					swap(&arr[j], &arr[j + 1]);
	}
#pragma omp single
	merge_sorted(arr, half, arr + half, size - half, dst);
}

int main(int argc, char* argv[]) {

	int progress = 0;
	const int A = 440;
	int N, num_threads;
	struct timeval T3, T4;
	double T1, T2, min, amount;


	if (argc < 3) {
		printf("Need to enter arguments (sizeArray) and num_threads");
		return -1;
	}

	N = atoi(argv[1]);
	num_threads = atoi(argv[2]);


#ifdef _OPENMP
	T1 = omp_get_wtime();
#else
	gettimeofday(&T3, NULL);
#endif
	
	unsigned int* restrict first_random = 0;
	unsigned int* restrict second_random = 0;

omp_set_nested(1);
#pragma omp parallel sections num_threads(2) shared(progress) 
{
#pragma omp section
	{

#ifdef _OPENMP
		omp_set_dynamic(0);
		omp_set_num_threads(num_threads);
#endif

		for (int i = 0; i < 100; i++) {

			double* restrict M1 = (double*)malloc(N * sizeof(double));
			double* restrict M2 = (double*)malloc(N / 2 * sizeof(double));
			double* restrict M2_Copy = (double*)malloc(N / 2 * sizeof(double));
			double* restrict dst = (double*)malloc(N / 2 * sizeof(double));

			for (int k = 0; k < N; k++) {
				srand(f(i));
				double value = 1 + rand_r(&first_random) % A;
				M1[k] = value;

			};


			for (int k = 0; k < N / 2; k++) {
				srand(f(i));
				double value = 1 + rand_r(&second_random) % (A * 10);
				M2[k] = value;
				M2_Copy[k] = value;

			}


#pragma omp parallel default(none) shared(N, M1, M2, M2_Copy, dst, min, amount)  
			{

				#pragma omp for 
				for (int k = 0; k < N; k++) {
					M1[k] = cosh(sqrt(M1[k])) / sinh(sqrt(M1[k]));
				}

				#pragma omp for  
				for (int k = 0; k < N / 2; k++) {
					if (k == 0) {
						M2[k] = fabsf(cos(M2[k]));
						continue;
					}
					M2[k] = fabsf(cos(M2[k] + M2_Copy[k - 1]));
				}
				#pragma omp barrier
				#pragma omp for
				for (int k = 0; k < N / 2; k++) {
					M2[k] = fabsf(M1[k] - M2[k]);
				}

				sort(M2, N / 2, dst);
		

				#pragma omp single
				min = dst[0];

				#pragma omp for reduction(+ : amount)
				for (int k = 0; k < N / 2; k++) {
					if (((long)(dst[k] / min) % 2) == 0) {
						amount += sin(dst[k]);
					}
				}
			}

			printf("Result: %f\n", amount);
			free(M1);
			free(M2);
			free(M2_Copy);
			free(dst);
			progress++;
		}

		double delta_ms;
#ifdef _OPENMP
		T2 = omp_get_wtime();
		delta_ms = T2 - T1;
#else
		gettimeofday(&T4, NULL);
		delta_ms = 1000 * (T4.tv_sec - T3.tv_sec) + (T4.tv_usec - T3.tv_usec) / 1000;
#endif
		printf("\nN=%d Milliseconds passed 100: %f\n", N, delta_ms);
	}

#pragma omp section 
	{
		double time = 0;
		while (progress < 100) {
			double time_temp = omp_get_wtime();
			if (time_temp - time < 1) {
				usleep(100);
				continue;
			};
			printf("%f%\n", (double)progress / 100 * 100.0);
			time = time_temp;
		}
	}
}
	return 0;
	
}