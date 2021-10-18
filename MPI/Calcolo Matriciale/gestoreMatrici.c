#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "gestoreMatrici.h"

int M =  3;
int N = 3;

int maxValue = 100;

void creaMatrice(double matrix[M][N]){
	srand(0);
	for(int i = 0; i < M ; i++){
		for(int j = 0; j < N; j++){	
			double value = (rand() > RAND_MAX / 2 ? -1 : 1) * ( (double)maxValue / RAND_MAX) * rand();
			matrix[i][j] = value;
			printf("%lf\n", matrix[i][j]);
			printf("i = %d, j = %d\n", i, j);
		}
	}
}

void scriviMatrice(double matrix[M][N])
{
	FILE *f = fopen("matrice", "w");
	
	int written = fwrite(matrix, sizeof(double), M*N, f);	
	if(written != M*N){
		perror("Scritti meno elementi di quanti se ne dovevano scrivere");
	}
	fflush(f);
	fclose(f);
	
	
}

void leggiMatrice(double matrix[M][N])
{
	FILE *f = fopen("matrice", "r");
	
	int read = fread(matrix, sizeof(double), M*N, f);	
	if(read != M*N){
		perror("Letto un numero di elementi di quanti se ne dovevano leggere");
	}
	fclose(f);
	
}

void leggiSottoMatrice(double *partialM, int startRow, int numRows)
{
	FILE *f = fopen("matrice", "r");
	int r = fseek(f, sizeof(double)*startRow*N, SEEK_SET);
	
	int read = fread(partialM, sizeof(double), numRows*N, f);	
	if(read != numRows*N){
		perror("Letto un numero di elementi diverso a quanti se ne dovevano leggere");
	}
	fclose(f);
	
}

void stampaMatrice(double matrix[M][N]){
	for(int i = 0; i < M ; i++){
		for(int j = 0; j < N; j++){
			printf("%lf\n", matrix[i][j]);
		}
	}
}

void scriviMatriceComeTxt(double matrix[M][N]){
	FILE *f = fopen("matrice.txt", "w");
	
	for(int i = 0; i < M ; i++){
		char *buffer = calloc(1000000, sizeof(char));
		int index = 0;
		for(int j = 0; j < N; j++){
			if(j < N -1){
				index += sprintf(buffer+index, "%lf ", matrix[i][j]);
			}else {
				index += sprintf(buffer+index, "%lf", matrix[i][j]);
			}
		}
		index += sprintf(buffer+index, "\n");
		printf("%ld\n%s\n", strlen(buffer), buffer);
		int written = fwrite(buffer, sizeof(char), strlen(buffer), f);	
		if(written != strlen(buffer)){
			perror("Scritti meno elementi di quanti se ne dovevano scrivere");
		}
		fflush(f);
	}
	fclose(f);
}
/*
int main(int argc, char* argv[])
{
	
	double matrix[M][N];
	
	creaMatrice(matrix);
	scriviMatrice(matrix);
	
	
	//leggiMatrice(matrix);
	//stampaMatrice(matrix);
	scriviMatriceComeTxt(matrix);
	
}
*/

void primaColonna(double colonna[M]){
	double matrix[M][N];
	leggiMatrice(matrix);
	
	for(int i = 0; i < M; i++){
		colonna[i] = matrix[i][0];
	}
}
