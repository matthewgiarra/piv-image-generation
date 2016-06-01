// Include the standard libraries.
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>


int main(int argc, char *argv[]){
	
	int thread_id;
	
#pragma omp parallel num_threads(5)
	{
		thread_id = omp_get_thread_num();
		
		printf("Hello from thread %d\n", thread_id);
	}
	
	return(0);
}
