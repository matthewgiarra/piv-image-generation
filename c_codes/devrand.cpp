/* This file demonstrates how to create random numbers in C using the random number generator file /dev/urandom. */

// Include the standard libraries.
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "velocity_functions.h"
#include "image_writing.h"
#include "tiffio.h"

// #include "myfunctions.h"

void devrand(double *array, int array_length, double lower_bound, double upper_bound) {

	// Initialize some variables
	FILE *urandom; // This will hold the file /dev/urandom
	unsigned int seed; // Random number seed
	
	int k; // Counter for loops

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
	
}

// Main function
int main(int argc, char * argv[]){
	
	// Function prototype
	// void plane_poiseuille(double *x, double *y, double *x_new, double h, double num_points);
	//
	// Array size
	int array_length = atoi(argv[1]);
	
	// Image dimensions
	const int image_height = 1024;
	const int image_width = 1024;
	
	// Coordinate system origin
	const double yc = image_height / 2;
	const double zc = 0;
	const double h = image_height / 2;
	
	// Max velocity 
	const double v_max = 20;
	
	//Allocate the array
	double *x 	  = (double*) malloc(array_length * sizeof(double));
	double *y 	  = (double*) malloc(array_length * sizeof(double));
	double *z 	  = (double*) malloc(array_length * sizeof(double));
	double *x_new = (double*) malloc(array_length * sizeof(double));
	
	// Error handling
	if(!x){
		std::cout << "Error: failed to allocate array x.\n";
		return(-1);
	}
	
	// Error handling
	if(!y){
		std::cout << "Error: failed to allocate array y.\n";
		return(-1);
	}
	
	// Error handling
	if(!z){
		std::cout << "Error: failed to allocate array z.\n";
		return(-1);
	}
	
	// Error handling
	if(!x_new){
		std::cout << "Error: failed to allocate array x_new.\n";
		return(-1);
	}
	
	// Pass array to the random number generator
	devrand(x, array_length, 0, image_width - 1);
	devrand(y, array_length, 0, image_height - 1);
	devrand(z, array_length, -1, 1);
			
	// Advect positions
	poiseuille(x, y, z, x_new, yc, zc, h, v_max, array_length);
	
	// Free the array
	free(x);
	free(y);
	free(z);
	free(x_new);
	
	// Create the file path to the saved image.
	std::string output_file_path_string_01 = "/Users/matthewgiarra/Desktop/image_01.tiff";
	// Create the file path to the saved image.
	std::string output_file_path_string_02 = "/Users/matthewgiarra/Desktop/image_02.tiff";
	
	// Specify a file path
	char *output_file_path_01 = (char*)output_file_path_string_01.c_str();
	char *output_file_path_02 = (char*)output_file_path_string_02.c_str();

	// Allocate memory for a 16-bit "slice" which will hold the image data.
	uint16_t *slice = new uint16_t[image_height * image_width];
	
	// Allocate memory for a linear index.
	int *linear_ind = (int*) malloc(sizeof(int));
	
	for(int k = 0; k < array_length; k ++){
		if((int)x[k] < image_width){
			// Calculate a linear index
			sub2ind(image_width, (int)y[k], (int)x[k], linear_ind);		
		
			// Set the pixel to white 
			slice[*linear_ind] = pow(2, 16) - 1;
		}
	}
		
	// Write a tiff image.
	writeTiff_bw16(output_file_path_01, slice, image_height, image_width);
	
	// Zero the slice
	for(int k = 0; k < image_height * image_width - 1; k ++){
		slice[k] = 0;
	}
	
	for(int k = 0; k < array_length; k ++){
		
		if((int)x_new[k] < image_width){
			// Calculate a linear index
			sub2ind(image_width, (int)y[k], (int)x_new[k], linear_ind);		
		
			// Set the pixel to white 
			slice[*linear_ind] = pow(2, 16) - 1;
		}
		
	}
		
	// Write a tiff image.
	writeTiff_bw16(output_file_path_02, slice, image_height, image_width);
	
	// GTFO
	return(0);
}

























