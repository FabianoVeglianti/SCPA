#include <mpi.h>
#include <stdio.h>

float prodottoPunto(float *x, float*y, int size){
	float partialResult = 0.0;
	for(int i = 0; i < size; i++){
		partialResult += *(x+i) * *(y+i);
	}
	return partialResult;
	
}

int main(int argc, char* argv[]){
	
	int p;
	int my_rank;
	int tag = 0;
	float result;
	
	MPI_Status status;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	
	float x[my_rank + 5], y[my_rank + 5];
	
	for(int i = 0; i < my_rank + 5; i++){
		x[i] = (float) i;
		y[i] = (float) i;
		/*
		if(i%2 == 0){
			y[i] = i;
		} else {
			y[i] = (float) -i;
		}
		*/
	}
	
	float partialResult = prodottoPunto(x,y, my_rank+5);
	
	MPI_Reduce(&partialResult, &result, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

	if(my_rank == 0){
		printf("Risultato = %f", result);
	}
	MPI_Finalize();
	return 0;
	
	
}
