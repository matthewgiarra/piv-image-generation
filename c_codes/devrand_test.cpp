/* This file demonstrates how to create random numbers in C using the random number generator file /dev/urandom. */

// Include the standard libraries.
#include <iostream>
#include <stdlib.h>
#include <math.h>
// #include "myfunctions.h"

int main(int argc, char * argv[]) {

	// Initialize some variables
	FILE *urandom; // This will hold the file /dev/urandom
	unsigned int seed; // Random number seed
	
	int k; // Counter for loops
	int lower_bound = 0;
	int upper_bound = 100;
	int number_of_numbers = 10;

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
	
	// Initialize x;
	double x[number_of_numbers];
	
	// Generate the numbers between the upper and lower bounds.
	for (k = 1; k < number_of_numbers; k++) {
		x[k] = rand() / (double)RAND_MAX * (upper_bound - lower_bound) + lower_bound ;
	}
	
	// Print the numbers
	for ( k = 1; k < number_of_numbers; k++ ){
		printf("%0.6f\n", x[k]);
	}

	printf("\n\n");

	// Return the pointer to the first element of the 
	// array containint the variables
	return(0);

}





