#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define maxInt 9
#define minInt 0

double scalarProduct(int *x, int *y, int size){
	double partialResult = 0.0;
	for(int i = 0; i < size; i++){
		partialResult += (double) x[i] * y[i];
	}
	return partialResult;
	
}

void computePartialY(double *partialY, int *partialM, int numRows, int *x, int N){
		for(int i = 0; i < numRows; i++){
			double value = scalarProduct((partialM+i*N), x, N);
			*(partialY+i) = value;			 
		}
}

void generateRandomIntVector(int *x, int M){
	srand(8555);
	for(int i = 0; i < M ; i++){
		x[i]  = (rand() % (maxInt - minInt + 1) + minInt);	
	}
}

void generateRandomIntMatrix(int *m, int M, int N, int my_rank){
	srand(my_rank*100);
	for(int i = 0; i < M ; i++){
		for(int j = 0; j < N; j++){
			*(m+i*N+j) = (rand() % (maxInt - minInt + 1) + minInt);	
		}
	}
}

int main(int argc, char* argv[])
{
	
	int p;
	int my_rank;
	int tag = 0;
	int M, N;
	
	MPI_Status status;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p); 

	int *x;
	double *y;

	/*Leggo le dimensioni della matrice*/
	if(my_rank == 0){
		printf("Inserisci #rows e #cols\n");
		scanf("%d %d", &M, &N);
		printf("\n");
		y = malloc(M*sizeof(double));
	}

	MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

	/*Il processo 0 genera il vettore x e ne fa il bcast*/
	x = malloc(N * sizeof(int));
	if(my_rank == 0){
		generateRandomIntVector(x, N);
		printf("x = ");
		for(int k = 0; k < N; k++){
			printf("%d ", x[k]);
		}
		printf("\n\n");
	}
	MPI_Bcast(x, N, MPI_INT, 0, MPI_COMM_WORLD);

	/*Ciascun processo calcola quali e quante righe della matrice deve gestire*/
	int r = M%p;
	int numRows, startRowIndex;
	if(my_rank < r){
		numRows = M/p + 1;
		startRowIndex = my_rank * numRows;
	} else {
		numRows = M/p +1;
		startRowIndex = r * numRows;
		numRows = M/p;
		startRowIndex += (my_rank - p + r) * numRows;
	}
//	printf("Process %d, startRowIndex = %d, numRows = %d\n", my_rank, startRowIndex, numRows);
	
	/* Per ogni processo la sottomatrice gestita corrisponde ad un certo numero di righe contigue della matrice */
	int partialM[numRows][N];
	generateRandomIntMatrix(&partialM[0][0], numRows, N, my_rank);
	
	/*Stampo la sottomatrice - per la validazione del codice*/
	char matrixString[1000];
	char *startString = matrixString;
	int w = sprintf(startString, "Processo %d\n", my_rank);
	startString = startString + w;
	for(int p = 0; p < numRows; p++){
		for(int q = 0; q < N; q++){
			w = sprintf(startString, "%d ", partialM[p][q]); 
			startString = startString + w;
		}
		w = sprintf(startString, "\n");
		startString = startString + w;
	}
	printf("%s\n", matrixString);

	/*Ogni processo calcola la parte del vettore y che gli compete*/
	double partialY[numRows];
	computePartialY(partialY, &partialM[0][0], numRows, x, N);
	
	/*
	for(int l = 0; l < numRows; l++){
		printf("Processo = %d, partialY[%d] = %lf\n", my_rank, l, partialY[l]);
	}
	*/
	
	/*Il processo 0 calcola quanti elementi di y si aspetta di ricevere da ciascun processo e dunque il displacement degli elementi ricevuti all'interno del vettore*/
	int elements = 0;
	int displs[p];
	int recvCounts[p];
	if(my_rank == 0){
		for(int k = 0; k < p; k ++){
			if(k < r){
				recvCounts[k] = M/p + 1;
				if(k == 0){
					displs[k] = 0;
				} else {
					displs[k] = displs[k-1] + recvCounts[k-1];
				}
			} else {
				recvCounts[k] = M/p;
				displs[k] = displs[k-1] + recvCounts[k-1];
			}
//			printf("recv[%d] = %d, displs[%d] = %d\n", k, recvCounts[k], k, displs[k]);
		}
	}
	
	
	/*Il processo 0 riceve i risultati parziali di y e li unisce in un unico vettore*/
	MPI_Gatherv(partialY, numRows, MPI_DOUBLE, y, recvCounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	/*Per fare in modo che l'ultima cosa stampata sia y*/
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);
	
	/*Il processo 0 stampa y - per la validazione*/
	if(my_rank == 0){
		printf("y = \n"); 
		for(int m = 0; m < N; m++){
			printf("%lf\n", *(y+m));
		}
		printf("\n");
	}
	
	MPI_Finalize();
	return 0;
}
