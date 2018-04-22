// classic.cpp : "Textbook" implementation of matrix multiply

// Author:  Paul J. Drongowski
// Address: Boston Design Center
//          Advanced Micro Devices, Inc.
//          Boxborough, MA 01719
// Date:    20 October 2005
//
// Copyright (c) 2005 Advanced Micro Devices, Inc.

// Celem tego programu jest prezentacja pomiaru i analizy 
//efektywnosci programu za pomoc¹  CodeAnalyst(tm).
// Implementacja mno¿enia macierzy jest realizowana za pomoca typowego 
// algorytmu podrêcznikowego. 
//  Efektywnoœæ tego podejscia jest niska poprzez 
// nieefektywn¹  kolejnosæ odwo³añ do elementów macierzy.
#include "stdafx.h"

#include <stdio.h>
#include <iostream>
#include <time.h>
#include <windows.h>
#include <math.h>
#include "omp.h"

#define MAXTHREADS 128
int NumThreads;
double start;

static const int ROWS = 3000;     // liczba wierszy macierzy
static const int COLUMNS = 3000;  // lizba kolumn macierzy

float matrix_a[ROWS][COLUMNS];    // lewy operand 
float matrix_b[ROWS][COLUMNS];    // prawy operand
float matrix_r[ROWS][COLUMNS];    // wynik

FILE *result_file;

void initialize_matrices()
{
	// zdefiniowanie zawarosci poczatkowej macierzy
#pragma omp parallel for 
	for (int i = 0; i < ROWS; i++) {
		for (int j = 0; j < COLUMNS; j++) {
			matrix_a[i][j] = (float) rand() / RAND_MAX;
			matrix_b[i][j] = (float) rand() / RAND_MAX;
			matrix_r[i][j] = 0.0;
		}
	}
}

void initialize_matricesZ()
{
	// zdefiniowanie zawarosci poczatkowej macierzy
#pragma omp parallel for 
	for (int i = 0; i < ROWS; i++) {
		for (int j = 0; j < COLUMNS; j++) {
			matrix_r[i][j] = 0.0;
		}
	}
}
void print_result()
{
	// wydruk wyniku
	for (int i = 0; i < ROWS; i++) {
		for (int j = 0; j < COLUMNS; j++) {
			fprintf(result_file, "%6.4f ", matrix_r[i][j]);
		}
		fprintf(result_file, "\n");
	}
}

void print_a()
{
	// wydruk wyniku
	for (int i = 0; i < ROWS; i++) {
		for (int j = 0; j < COLUMNS; j++) {
			fprintf(result_file, "%6.4f ", matrix_a[i][j]);
		}
		fprintf(result_file, "\n");
	}
}

void print_b()
{
	// wydruk wyniku
	for (int i = 0; i < ROWS; i++) {
		for (int j = 0; j < COLUMNS; j++) {
			fprintf(result_file, "%6.4f ", matrix_b[i][j]);
		}
		fprintf(result_file, "\n");
	}
}

void multiply_matrices_IJK_IKJ(int r)
{
	// mnozenie macierzy 
	for (int i = 0; i < ROWS; i+=r) {
		for (int j = 0; j < COLUMNS; j+=r) {
			for (int k = 0; k < COLUMNS; k+=r) {
#pragma omp parallel for 
				for (int ii = i; ii < i + r; ii++) {
					for (int kk = k; kk < k + r; kk++) {
						for (int jj = j; jj < j + r; jj++) {
							matrix_r[ii][jj] += matrix_a[ii][kk] * matrix_b[kk][jj];
						}
					}
				}
			}
		}
	}
}

void multiply_matrices_JKI()
{
	// mnozenie macierzy 
#pragma omp parallel for 
	for (int j = 0; j < COLUMNS; j++)
		for (int k = 0; k < COLUMNS; k++)
			for (int i = 0; i < ROWS; i++)
				matrix_r[i][j] += matrix_a[i][k] * matrix_b[k][j];

}


void print_elapsed_time()
{
	double elapsed;
	double resolution;

	// wyznaczenie i zapisanie czasu przetwarzania
	elapsed = (double)clock() / CLK_TCK;
	resolution = 1.0 / CLK_TCK;

	printf("Czas wykonania programu: %8.4f sec (%6.4f sec rozdzielczosc pomiaru)\n",
		elapsed - start, resolution);
}

int main(int argc, char* argv[])
{
	//	 start = (double) clock() / CLK_TCK ;

	bool USE_MULTIPLE_THREADS = false;

	if (argc > 1) {
		USE_MULTIPLE_THREADS = true;
	}

	//Determine the number of threads to use
	if (USE_MULTIPLE_THREADS) {
		SYSTEM_INFO SysInfo;
		GetSystemInfo(&SysInfo);
		NumThreads = SysInfo.dwNumberOfProcessors;
		if (NumThreads > MAXTHREADS)
			NumThreads = MAXTHREADS;
	}
	else
		NumThreads = 1;
	printf("Klasyczny algorytm mnozenia macierzy, liczba watkow %d \n", NumThreads);
	printf("liczba watkow  = %d\n\n", NumThreads);

	int r = (int)sqrt((double)6*1024*1024 / sizeof(float) / 3);
	r = ROWS / ceil((float) ROWS / r);

	omp_set_num_threads(NumThreads);

	initialize_matrices();
	start = (double)clock() / CLK_TCK;
	multiply_matrices_IJK_IKJ(r);
	print_elapsed_time();

	/*result_file = fopen("abc.txt", "w");

	print_a();
	print_b();
	print_result();*/
	

	initialize_matricesZ();
	start = (double)clock() / CLK_TCK;
	multiply_matrices_JKI();
	print_elapsed_time();

	/*fclose(result_file);*/
	//system("pause");
	return(0);
}