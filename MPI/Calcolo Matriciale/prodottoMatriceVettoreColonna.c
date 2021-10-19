#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define maxInt 9
#define minInt 0

double scalarProduct(int *x, int *y, int size, int startIndex){
	double partialResult = 0.0;
	for(int i = 0; i < size; i++){
		partialResult += (double) x[i] * y[startIndex + i];
	}
	return partialResult;
	
}

void computeContributeToY(double *contibuteToY, int *partialM, int M, int *x,int N,  int numCols, int startIndex){
		for(int i = 0; i < M; i++){
			double value = scalarProduct((partialM + i*numCols), x, numCols, startIndex);
			*(contibuteToY+i) = value;			 
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

	/*Leggo la dimensione della matrice*/
	if(my_rank == 0){
		printf("Inserisci #rows e #cols\n");
		scanf("%d %d", &M, &N);
		printf("\n");
		y = malloc(M*sizeof(double));
	}

	MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

	/*Genero x e ne faccio il bcast*/
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

	/*Ciascun processo si calcola quali e quante colonne della matrice generare*/
	int r = N%p;
	int numCols, startColsIndex;
	if(my_rank < r){
		numCols = N/p + 1;
		startColsIndex = my_rank * numCols;
	} else {
		numCols = N/p +1;
		startColsIndex = r * numCols;
		numCols = N/p;
		startColsIndex += (my_rank - p + r) * numCols;
	}
//	printf("Process %d, startRowIndex = %d, numRows = %d\n", my_rank, startRowIndex, numRows);
	
	/* Per processo la sottomatrice gestita corrisponde ad un certo numero di colonne contigue della matrice */
	int partialM[M][numCols];
	generateRandomIntMatrix(&partialM[0][0], M, numCols, my_rank);
	
	/*Stampo la sottomatrice - per la validazione del codice*/
	char matrixString[1000];
	char *startString = matrixString;
	int w = sprintf(startString, "Processo %d\n", my_rank);
	startString = startString + w;
	for(int p = 0; p < M; p++){
		for(int q = 0; q < numCols; q++){
			w = sprintf(startString, "%d ", partialM[p][q]); 
			startString = startString + w;
		}
		w = sprintf(startString, "\n");
		startString = startString + w;
	}
	printf("%s\n", matrixString);

	/*Calcolo il contributo a y del processo*/
	double contributeToY[M];
	computeContributeToY(contributeToY, &partialM[0][0], M, x, N, numCols, startColsIndex);
	
	/*Stampo il contributo a y del processo - per la validazione*/
	for(int l = 0; l < M; l++){
		printf("Processo = %d, partialY[%d] = %lf\n", my_rank, l, contributeToY[l]);
	}
		
	/*Sommo i contributi a y di ciascun processo*/
	MPI_Reduce(contributeToY, y, M,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	
	/*Per fare in modo che l'ultima cosa stampata sia y, MA NON FUNZIONA*/
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
