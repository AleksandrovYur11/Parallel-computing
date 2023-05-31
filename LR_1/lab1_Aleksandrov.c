#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>

void func_map_m1(int sizeArrayM1, double M1[]) {
	for (int k = 0; k < sizeArrayM1; k++) {
		M1[k] = cosh(sqrt(M1[k]))/sinh(sqrt(M1[k]));
	}
}
 
void func_map_m2(int sizeArrayM2, double M2[]) {
	double cloneArrayM2[sizeArrayM2];
	for (int k = 0; k < sizeArrayM2; k++) {
		cloneArrayM2[k] = M2[k];
	}

	for (int k = 0; k < sizeArrayM2; k++) {
		if (k == 0) {
			M2[k] = fabsf(cos(M2[k]));
			continue;
		}
		M2[k] = fabsf(cos(M2[k] + cloneArrayM2[k - 1]));
	}
}

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

double reduce(double M2[], int sizeArrayM2) {
	double amount = 0;
	double min = 0;
	for (int i = 0; i < sizeArrayM2; i++) {
		if (M2[i] != 0) {
			min = M2[i];
			break;
		}
	}

	for (int k = 0; k < sizeArrayM2; k++) {
		if (((long)(M2[k] / min) % 2) == 0) {
			amount += sin(M2[k]);
		}
	}
	return amount;
}


int main(int argc, char* argv[]) {

	const int A = 440;
	int N;
	struct timeval T1, T2;
	long delta_ms;

	if (argc < 2) {
		printf("Need to enter arguments (sizeArray)");
		return -1;
	}

	N = atoi(argv[1]);
	
	gettimeofday(&T1, NULL);

	unsigned int first_random = 0;
	unsigned int second_random = 0;

	for (int i = 0; i < 100; i++) {

		srand(i);

		double* M1 = (double*)malloc(N * sizeof(double));
		double* M2 = (double*)malloc(N / 2 * sizeof(double));

		for (int k = 0; k < N; k++) {
			double value = 1 + rand_r(&first_random) % A;
			M1[k] = value;
			//printf("%.2f\n", value);
		};

		for (int k = 0; k < N / 2; k++) {
			double value = 1 + rand_r(&second_random) % (A * 10);
			M2[k] = value;
			//printf("%.2f\n", value);
		}

		func_map_m1(N, M1);

		func_map_m2(N / 2, M2);
		
		merge(M1, M2, N / 2);

		//stupid_sort(M2, N / 2);

		int result = reduce(M2, N / 2);

		printf("Result: %d\n", result);

		free(M1);
		free(M2);
	}


	gettimeofday(&T2, NULL);
	delta_ms = 1000 * (T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec) / 1000;
	printf("\nN=%d. Milliseconds passed: %ld\n", N, delta_ms);

	return 0;
}
