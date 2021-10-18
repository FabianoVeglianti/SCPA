#include <stdio.h>
#include <string.h>
#include <mpi.h>

float f(float x){
	return x;
}
	
float Trap(float local_a, float local_b, int local_n, float h)
{
	float x, integral;
	int i;
	
	integral = (f(local_a)+f(local_b))/2.0;
	for(i = 1; i < local_n; i++){
		x = local_a + i*h;
		integral +=f(x);
	}
	integral = integral * h;
	return integral;
}



int main(int argc, char* argv[])
{
	int my_rank;
	int p;
	int i;
	int tag = 0;
	float a,b;
	int n;
	MPI_Status status;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	
	if(my_rank == 0){
		printf("Inserisci a, b e n\n");
		scanf("%f %f %d", &a, &b, &n);
	}
	MPI_Bcast(&a, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&b, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	float h, local_a, local_b, local_int, integral;;
	int local_n;
	int r;
	
	h=(b-a)/n;
	r = n%p;
	
	/*
	
	if(my_rank >= p-r && my_rank < p-1){
		local_n = n/p + 1;
		local_a = a + (p-1-r)*(local_n -1)*h + (my_rank - p + r)*local_n*h;
	} else {
		local_n = n/p;
		local_a = a + my_rank * local_n * h;
	}
	
	local_b = local_a + local_n * h;
	
	*/
	
	if(my_rank < r){
		local_n = n/p + 1;
		local_a = a + my_rank * local_n * h;
	} else {
		local_n = n/p + 1; local_a = a + r * local_n * h;
		local_n = n/p; local_a = local_a + (my_rank - r)*local_n * h;
	}
	local_b = local_a + local_n * h;
	
	
	local_int = Trap(local_a, local_b, local_n, h);
	
	MPI_Reduce(&local_int, &integral, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	if(my_rank == 0){
		printf("Con n=%d trapezzi, stimiamo l'integrale da %f a %f pari a %f",n, a,b,integral);
	}
	MPI_Finalize();
	return 0;
	
}
