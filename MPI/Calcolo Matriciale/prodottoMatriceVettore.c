#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define maxInt 9
#define minInt 0

double scalarProduct(double *x, double *y, int size){
	double partialResult = 0.0;
	for(int i = 0; i < size; i++){
		partialResult += (double) x[i] * y[i];
	}
	return partialResult;
	
}

void computeContributeToPartialY(double *contributeToPartialY, double *partialM, int numRows, int numCols, double *x){
		for(int i = 0; i < numRows; i++){
			double value = scalarProduct((partialM+i*numCols), x, numCols);
			*(contributeToPartialY+i) = value;			 
		}
}

void generateRandomIntVector(double *x, int M, int seed){
	srand(seed);
	for(int i = 0; i < M ; i++){
		x[i]  = (double)(rand() % (maxInt - minInt + 1) + minInt);	
	}
}

void generateRandomIntMatrix(double *m, int M, int N, int my_rank){
	srand(my_rank*100);
	for(int i = 0; i < M ; i++){
		for(int j = 0; j < N; j++){
			*(m+i*N+j) = (double)(rand() % (maxInt - minInt + 1) + minInt);	
		}
	}
}

int main(int argc, char* argv[])
{
	
	int p;
	int my_rank;
	int tag = 0;
	int M, N, iterations;
	
	MPI_Status status;
	
	MPI_Init(&argc, &argv);
	
	MPI_Comm newComm;
	MPI_Comm_dup(MPI_COMM_WORLD, &newComm);
	
	MPI_Comm_rank(newComm, &my_rank);
	MPI_Comm_size(newComm, &p); 

	double *x;
	double *y;

	if(my_rank == 0){
		printf("Inserisci #rows, #cols e #iterations\n");
		scanf("%d %d %d", &M, &N, &iterations);
		printf("\n");
		y = calloc(M, sizeof(double));
	}

	MPI_Bcast(&M, 1, MPI_INT, 0, newComm);
	MPI_Bcast(&N, 1, MPI_INT, 0, newComm);
	MPI_Bcast(&iterations, 1, MPI_INT, 0, newComm);

	/*Stabilisco quali colonne un processo ha in base al rank*/
	int row_r = M%2;
	int numRows, startRowIndex;
	int orizzontalId;
	if(my_rank < p/2){
		orizzontalId = 0;
	} else {
		orizzontalId = 1;
	}
	if(orizzontalId < row_r){
		numRows = M/2 + 1;
		startRowIndex = orizzontalId * numRows;
	} else {
		numRows = M/2 +1;
		startRowIndex = row_r * numRows;
		numRows = M/2;
		startRowIndex += (orizzontalId - row_r) * numRows;
	}
	
	/*Stabilisco quali colonne un processo ha in base al rank*/
	int processPerRow = p/2; // ogni riga è divisa tra p/2 processi
	int col_r = N%processPerRow;
	int numCols, startColsIndex;
	int verticalId=my_rank%processPerRow;
	if(verticalId < col_r){
		numCols = N/processPerRow + 1;
		startColsIndex = verticalId * numCols;
	} else {
		numCols = N/processPerRow +1;
		startColsIndex = col_r * numCols;
		numCols = N/processPerRow;
		startColsIndex += (verticalId - col_r) * numCols;
	}

	
	double *contributeToPartialY;
	double *partialY;

	/*Per ogni gruppo di colonne definisco un gruppo e un comunicatore
	 * Ad esempio una matrice 8x8 con p=6 in cui la matrice è divisa come segue
	 * 	0		1		2
	 * 	3		4		5
	 * i gruppi di colonne sono (0,3); (1,4) e (2,5)
	 * Ciò è fatto in base al resto della divisione tra il rango del processo e il numero di processi in una riga.
	 * 
	 * */
	MPI_Group all;
	MPI_Comm_group(newComm, &all);
	MPI_Group colsGroup;
	MPI_Comm colsComm;
	
	int colRanks[2];
	int groupId = my_rank%processPerRow;
	int counter = 0;
	for(int rank = 0; rank < p; rank++){
		if(rank%processPerRow == groupId){
			colRanks[counter] = rank;
			counter = counter +1;
		}
	}
	printf("Processo %d - %d %d\n", my_rank, colRanks[0], colRanks[1]);
	
	MPI_Group_incl(all, 2, colRanks, &colsGroup);
	MPI_Comm_create(newComm, colsGroup, &colsComm);
	
	
	/*Per ogni gruppo di colonne definisco un gruppo e un comunicatore
	 * Ad esempio una matrice 8x8 con p=6 in cui la matrice è divisa come segue
	 * 	0		1		2
	 * 	3		4		5
	 * i gruppi di colonne sono (0,1,2) e (3,4,5)
	 * Ciò è fatto (avendo fissato 2 come numero di gruppi) vedendo se il rango del processo è minore a p/2 oppure no
	 * 
	 * */
	counter = 0;
	int rowRanks[p/2];
	groupId = my_rank<p/2;
	for(int processRank = 0; processRank < p; processRank++){
		if((processRank<p/2)==groupId){
			rowRanks[counter] = processRank;
			counter = counter+1;
		}
	}

	MPI_Group rowsGroup;
	MPI_Comm rowsComm;
	MPI_Group_incl(all, p/2, rowRanks, &rowsGroup);
	MPI_Comm_create(newComm, rowsGroup, &rowsComm);
	
	
	int my_rank_col;
	int p_col;
	
	MPI_Comm_rank(colsComm, &my_rank_col);
	MPI_Comm_size(colsComm, &p_col); 
	
	int my_rank_row;
	int p_row;
	
	MPI_Comm_rank(rowsComm, &my_rank_row);
	MPI_Comm_size(rowsComm, &p_row); 
	
	printf("Processo %d - Processo in col %d - totale %d - Processo in row %d - totale %d\n", my_rank, my_rank_col, p_col, my_rank_row, p_row);

	printf("Processo %d\nStartRow %d - NumRows %d\nStartCol %d - NumCols %d\n", my_rank, startRowIndex, numRows, startColsIndex, numCols);
	
	
	for(int j = 0;  j < iterations; j++){
		
		if(j ==0) {
			x = malloc(numCols * sizeof(double));
			if(my_rank_col == 0){
				generateRandomIntVector(x, numCols, my_rank*78);
			}
			MPI_Bcast(x, numCols, MPI_DOUBLE, 0, colsComm);
			
			contributeToPartialY = malloc(numRows * sizeof(double));
			partialY = malloc(numRows * sizeof(double));
		} else {
			MPI_Reduce(contributeToPartialY, partialY, numRows, MPI_DOUBLE, MPI_SUM, 0, rowsComm);
			
			/*
			char t[1000];
			char *pt = t;
			int bw = sprintf(pt, "X dopo reduce Processo %d: ", my_rank);
			pt = pt + bw;
			for(int d= 0; d < numRows; d++){
				bw = sprintf(pt, "%lf ", contributeToPartialY[d]);
				pt = pt + bw;
			}
			bw = sprintf(pt, "\n");
			pt = pt + bw;
			printf("%s\n", t);
			*/
			
			sleep(5);
			
			if(my_rank%processPerRow == 0){
	
				int elements = 0;
				int displsGather[2];
				int recvCountsGather[2];
		
					
				if(my_rank == 0){
					for(int k = 0; k < 2; k ++){
						if(k < row_r){
							recvCountsGather[k] = M/2+ 1;
							if(k == 0){
								displsGather[k] = 0;
							} else {
								displsGather[k] = displsGather[k-1] + recvCountsGather[k-1];
							}
						} else {
							recvCountsGather[k] = M/2;
							if(k == 0){
								displsGather[k] = 0;
							} else {
								displsGather[k] = displsGather[k-1] + recvCountsGather[k-1];
							}
						}
			//			printf("recv[%d] = %d, displs[%d] = %d\n", k, recvCounts[k], k, displs[k]);
					}
				}	
			
				MPI_Gatherv(partialY, numRows, MPI_DOUBLE, y, recvCountsGather, displsGather, MPI_DOUBLE, 0, colsComm);
				
				if(my_rank == 0){
						printf("Vettore x = [%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f]\n", y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7]);
				}
				
			}
			
			if(my_rank_col == 0){
				
				int elements = 0;
				int displsScatter[p/2];
				int sendCountsScatter[p/2];
		
					
				if(my_rank == 0){
					for(int k = 0; k < processPerRow; k ++){
						if(k < col_r){
							sendCountsScatter[k] = N/processPerRow+ 1;
							if(k == 0){
								displsScatter[k] = 0;
							} else {
								displsScatter[k] = displsScatter[k-1] + sendCountsScatter[k-1];
							}
						} else {
							sendCountsScatter[k] = N/processPerRow;
							if(k == 0){
								displsScatter[k] = 0;
							} else {
								displsScatter[k] = displsScatter[k-1] + sendCountsScatter[k-1];
							}
						}
			//			printf("recv[%d] = %d, displs[%d] = %d\n", k, recvCounts[k], k, displs[k]);
					}
					
					for(int y = 0; y < processPerRow; y++){
						printf("%d %d %d\n", y, sendCountsScatter[y], displsScatter[y]);
					}
				}	
				sleep(5);
				MPI_Scatterv(y, sendCountsScatter, displsScatter, MPI_DOUBLE, x, numCols, MPI_DOUBLE, 0, rowsComm);
				/*
				char t[1000];
				char *pt = t;
				int bw = sprintf(pt, "X dopo scatter Processo %d: ", my_rank);
				pt = pt + bw;
				for(int d= 0; d < numCols; d++){
					bw = sprintf(pt, "%lf ", x[d]);
					pt = pt + bw;
				}
				bw = sprintf(pt, "\n");
				pt = pt + bw;
				printf("%s\n", t);
				*/
				
			}
			MPI_Bcast(x, numCols, MPI_DOUBLE, 0, colsComm);
			
		}
		
		sleep(3);
		char vector[300];
		char *startVector = vector;
		int wr = sprintf(startVector, "Vector Processo %d\n", my_rank);
		startVector = startVector + wr;
		for(int o = 0; o < numCols; o++){
			wr = sprintf(startVector, "%lf ", x[o]);
			startVector = startVector + wr;
		}
		wr = sprintf(startVector, "\n");
		startVector = startVector + wr;
		printf("%s\n", vector);
	//	printf("Process %d, startRowIndex = %d, numRows = %d\n", my_rank, startRowIndex, numRows);
		sleep(3);
		

		//ogni processo genera la sua sottomatrice di dimensione numRows x numCols
		double partialM[numRows][numCols];
		generateRandomIntMatrix(&partialM[0][0], numRows, numCols, my_rank);
		
		
		char matrixString[1000];
		char *startString = matrixString;
		int w = sprintf(startString, "Processo %d\n", my_rank);
		startString = startString + w;
		for(int r = 0; r < numRows; r++){
			for(int q = 0; q < numCols; q++){
				w = sprintf(startString, "%lf ", partialM[r][q]); 
				startString = startString + w;
			}
			w = sprintf(startString, "\n");
			startString = startString + w;
		}
		printf("%s\n", matrixString);
		

		computeContributeToPartialY(contributeToPartialY, &partialM[0][0], numRows, numCols, x);
		
		
		/*
		char contribute[1000];
		char *startContribute = contribute;
		int num = sprintf(startContribute, "Contributo Processo %d: ", my_rank);
		startContribute = startContribute + num;
		
		for(int t = 0; t < numRows; t++){
			num = sprintf(startContribute, "%lf ", contributeToPartialY[t]); 
			startContribute = startContribute + num;
		}
		num = sprintf(startContribute, "\n"); 
		startContribute = startContribute + num;
		printf("%s\n", contribute);
		sleep(3);
		* */
	}
	sleep(3);
	printf("Processo %d - ha finito\n", my_rank);
	sleep(3);
	MPI_Barrier(newComm);
	return 0;
	/*
		for(int l = 0; l < numRows; l++){
			printf("Processo = %d, partialY[%d] = %lf\n", my_rank, l, partialY[l]);
		}
		
		*/
		

//	MPI_Gatherv(y, numRows, MPI_DOUBLE, y, recvCounts, displs, MPI_DOUBLE, 0, newComm);
	MPI_Barrier(newComm);
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
