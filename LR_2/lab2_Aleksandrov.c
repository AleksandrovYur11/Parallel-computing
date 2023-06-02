#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <time.h>
#include "FW_1.3.1_Lin64/fwBase.h"
#include "FW_1.3.1_Lin64/fwSignal.h"


void func_map_m1(int sizeArrayM1, float M1[]) {

	float* arrayCosh = (float*)malloc(sizeArrayM1 * sizeof(float));
	float* arraySinh = (float*)malloc(sizeArrayM1 * sizeof(float));

	fwsSqrt_32f(M1, M1, sizeArrayM1);
	fwsCosh_32f_A11(M1, arrayCosh, sizeArrayM1);
	fwsSinh_32f_A11(M1, arraySinh, sizeArrayM1);
	fwsDiv_32f(arrayCosh, arraySinh, M1, sizeArrayM1);

	free(arraySinh);
	free(arrayCosh);
}

void func_map_m2(int sizeArrayM2, float M2[], float M1[]) {

	float* M2_tmp = (float*)malloc(sizeArrayM2 * sizeof(float));

	fwsSqrt_32f_A11(M1, M1, sizeArrayM2 * 2);
	fwsCos_32f_A11(M1, M1, sizeArrayM2 * 2);
	fwsInv_32f_A11(M1, M1, sizeArrayM2 * 2);
	fwsAdd_32f(M2_tmp, &M2_tmp[1], &M2[1], sizeArrayM2 - 1);
	fwsCos_32f_A11(M2, M2, sizeArrayM2);
	fwsAbs_32f_I(M2, sizeArrayM2);

	free(M2_tmp);
}

void merge(float M1[], float M2[], int sizeArrayM2) {
	fwsSub_32f(M1, M2, M2, sizeArrayM2);
	fwsAbs_32f(M2, M2, sizeArrayM2);
}


int isSorted(float M2[], int sizeArrayM2)
{
	while (--sizeArrayM2 > 1)
		if (M2[sizeArrayM2] < M2[sizeArrayM2 - 1])
			return 0;
	return 1;
}

void swap(float* x, float* y)
{
	float tmp = *x;
	*x = *y;
	*y = tmp;
}

void shuffle(float M2[], int sizeArrayM2)
{
	for (int i = 0; i < sizeArrayM2; i++)
		swap(&M2[i], &M2[rand() % sizeArrayM2]);
}

void stupid_sort(float M2[], int sizeArrayM2)
{
	while (!isSorted(M2, sizeArrayM2))
		shuffle(M2, sizeArrayM2);
}


double reduce(float M2[], int sizeArrayM2) {
	float amount = 0;
	float min = 0;
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

	fwSetNumThreads(numberFlow);

	gettimeofday(&T1, NULL);

	unsigned int first_random = 0;
	unsigned int second_random = 0;




	for (int i = 0; i < 100; i++) {

		srand(i);

		float* M1 = (float*)malloc(N * sizeof(float));
		float* M2 = (float*)malloc(N / 2 * sizeof(float));


		for (int k = 0; k < N; k++) {
			float value = 1 + rand_r(&first_random) % A;
			M1[k] = value;
	
		};

		for (int k = 0; k < N / 2; k++) {
			float value = 1 + rand_r(&second_random) % (A * 10);
			M2[k] = value;
		}


		func_map_m1(N, M1);

		func_map_m2(N / 2, M2, M1);

		merge(M1, M2, N / 2);

		stupid_sort(M2, N / 2);

		float result = reduce(M2, N / 2);

		printf("Result: %f\n", result);

		free(M1);
		free(M2);
	}


	gettimeofday(&T2, NULL);
	delta_ms = 1000 * (T2.tv_sec - T1.tv_sec) + (T2.tv_usec - T1.tv_usec) / 1000;
	printf("\nN=%d. Milliseconds passed: %ld\n", N, delta_ms);

	return 0;
}