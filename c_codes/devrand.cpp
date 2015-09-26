/* This file demonstrates how to create random numbers in C using the random number generator file /dev/urandom. */

// Include the standard libraries.
#include <iostream>
#include <stdlib.h>
#include <math.h>
// #include "myfunctions.h"

void devrand(double *array, int array_length) {

	// Initialize some variables
	FILE *urandom; // This will hold the file /dev/urandom
	unsigned int seed; // Random number seed
	
	int k; // Counter for loops
	int lower_bound = 0;
	int upper_bound = 100;

	// Inform the user
	std::cout << "Lower bound = " << lower_bound << ", upper bound = " << upper_bound << "\n";	
		
	// Open the /urandom file for reading.
	urandom = fopen("/dev/urandom", "r");

	// Inform the user if the file wasn't opened successfully.
	if (urandom == NULL) {
		fprintf(stderr, "Cannot open /dev/urandom!\n");
		exit(1);
	}

	// Put a random number from /urandom into the memory location of the seed; 
	// i.e. set the seed to a random number obtained from /urandom.
	// "&seed" specifies a pointer to the memory location of the variable 'seed';
	// sizeof(seed) specifies the size in bytes of the the variable "seed"
	// The number "1" here means "read 1 element from the file" (which is /urandom)
	// The last input is the pointer to the file /urandom (specified as a pointer in the beginning of the file)
	fread(&seed, sizeof(seed), 1, urandom);

	// Seed the random number generator using the seed obtained from /urandom
	srand(seed); 
	
	// Generate the numbers between the upper and lower bounds.
	for (k = 0; k < array_length; k++) {
		array[k] = rand() / (double)RAND_MAX * (upper_bound - lower_bound) + lower_bound ;
	}
	
	// Print the numbers
	std::cout << "Last number in array (internal function): \n";
	printf("%06.6f\n", array[array_length-1]); 

}

// Main function
int main(int argc, char * argv[]){
	
	// Array size
	int array_length = atoi(argv[1]);
	
	//Allocate the array
	double *array = (double*) malloc(array_length * sizeof(double));
	
	// Error handling
	if(!array){
		std::cout << "Error: failed to allocate array.\n";
		return(-1);
	}
	
	// Pass array to the random number generator
	devrand(array, array_length);
	
	// Print the numbers
	std::cout << "\n\nReading last number in array (from memory):\n";
	printf("%06.6f\n", array[array_length-1]); 
	
	// Free the array
	free(array);
	
	// GTFO
	return(0);
	
}

























