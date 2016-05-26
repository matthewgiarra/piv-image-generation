/* This file demonstrates how to create random numbers in C using the random number generator file /dev/urandom. */

// Include the standard libraries.
#include <iostream>
#include <stdlib.h>
#include <math.h>
// #include "velocity_functions.h"
// #include "image_writing.h"
// #include "tiffio.h"

// #include "myfunctions.h"
void generateParticleImage(double *B, double *X, double *Y, double *dp, double *I, int M, int N, int num_particles) /* Input variables */
{
	
	// INPUTS
	
	// nlhs is the number of output (left-side) arguments, or the size of the plhs array.
	
	// plhs is the array of output arguments.
	
	// nrhs is the number of input (right hand side) arguments, or the size of the prhs array.
	
	// prhs is the array of input arguments.
	
	// Pointers to arrays
    const double pi = 3.141592653589793;
	const double sqrt8 = sqrt(8);
	
	// Height and width of the image to be rendered
	int p, r, c, ind; // counters
	
	// Numbers of rows and columns
	M = *num_rows;
	N = *num_cols;
	
	// Define the cutoff intensity
	double cutoff_intensity = exp(-3);
	
	// Fractional rows and columns
	//to be rendered for each particle
	double minRenderedCol_fract;
	double maxRenderedCol_fract;
	double minRenderedRow_fract;
	double maxRenderedRow_fract;
	
	// Integer rows and columns
	//to be rendered for each particle
	int minRenderedCol;
	int maxRenderedCol;
	int minRenderedRow;
	int maxRenderedRow;
				
	// Loop over all the particles
	// Loop over all the particles
	for(p = 0; p < num_particles; p++){		
		
		// Min and max rows and columns to render for that column.
		minRenderedCol_fract = fmax(0, X[p] - 1 * dp[p]);
		maxRenderedCol_fract = fmin(N-1, X[p] + 1 * dp[p]);
		minRenderedRow_fract = fmax(0, Y[p] - 1 * dp[p]);
		maxRenderedRow_fract = fmin(M-1, Y[p] + 1 * dp[p]);
		
		// Integer rows to render
		minRenderedRow = (int) (minRenderedRow_fract - 
			fmod(minRenderedRow_fract, 1));
		maxRenderedRow = (int) (maxRenderedRow_fract - 
			fmod(maxRenderedRow_fract, 1));
	
		// Integer columns to render
		minRenderedCol = (int) (minRenderedCol_fract - 
			fmod(minRenderedCol_fract, 1));
		maxRenderedCol = (int) (maxRenderedCol_fract - 
			fmod(maxRenderedCol_fract, 1));
		
		// Skip particles that are outside of the domain
		if(minRenderedRow >=0 & maxRenderedRow < M & minRenderedCol >=0 & maxRenderedCol < N & I[p] > cutoff_intensity){
	
			printf("R_min = %d\t R_max = %d\nC_min = %d\tC_max = %d\n", minRenderedRow, maxRenderedRow, minRenderedCol, maxRenderedCol);
			
			// Loop over all the particles.
			for(r = minRenderedRow; r <= maxRenderedRow; r++){
				for(c = minRenderedCol; c <= maxRenderedCol; c++){
					
					// Index of the pixel to render					
					ind = r + M * c;

					// Add the intensity to the image
					B[ind] += dp[p] * 
						I[p] * dp[p] * dp[p] * pi / 32 *
                   (erf( sqrt8 *  (c - X[p] + 0.5) 
                   / dp[p]) - erf(sqrt8 * 
                   (c - X[p] - 0.5) / dp[p])) * 
                   (erf( sqrt8 *  (r - Y[p] + 0.5) 
                   / dp[p]) - erf(sqrt8 * 
                   (r - Y[p] - 0.5) / dp[p]));					
				}
			}			
		}		
	}
	
    return;
}




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
	
	
	// Image dimensions
	const int image_height = atoi(argv[1]);
	const int image_width  = atoi(argv[2]);
	int num_particles 	   = atoi(argv[3]);
		
	//Allocate memory for the particle positions and sizes
	double *x 	  = (double*) malloc(num_particles * sizeof(double));
	double *y 	  = (double*) malloc(num_particles * sizeof(double));
	double *dp	  = (double*) malloc(num_particles * sizeof(double));
	double *I	  = (double*) malloc(num_particles * sizeof(double));
	
	// Allocate memory for the image
	double *B = (double*) malloc(image_width * image_height * sizeof(double));
	
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
	if(!B){
		std::cout << "Error: failed to allocate array B.\n";
		return(-1);
	}
	
	// Error handling
	if(!dp){
		std::cout << "Error: failed to allocate array dp.\n";
		return(-1);
	}
	
	// Error handling
	if(!I){
		std::cout << "Error: failed to allocate array I\n";
		return(-1);
	}
	
	// Pass array to the random number generator
	// to generate particle positions
	devrand(x, num_particles, 0, image_width - 1);
	devrand(y, num_particles, 0, image_height - 1);
	
	// Random intensities and diameters
	devrand(I, num_particles, 0, 1);
	devrand(dp, num_particles, 1.5, 4.0);
	//
	// for(int p = 0; p < num_particles; p++){
	// 	I[p] = 1.00;
	// 	dp[p] = 2.80;
	// }
		
	// // Create the file path to the saved image.
	// std::string output_file_path_string_01 = "/Users/matthewgiarra/Desktop/image_01.tiff";
	// // Create the file path to the saved image.
	// std::string output_file_path_string_02 = "/Users/matthewgiarra/Desktop/image_02.tiff";
	//
	// // Specify a file path
	// char *output_file_path_01 = (char*)output_file_path_string_01.c_str();
	// char *output_file_path_02 = (char*)output_file_path_string_02.c_str();
	generateParticleImage(B, x, y, dp, I, image_height, image_width, num_particles);

	// // Allocate memory for a 16-bit "slice" which will hold the image data.
// 	uint16_t *slice = new uint16_t[image_height * image_width];
//
// 	// Write a tiff image.
// 	writeTiff_bw16(output_file_path_01, slice, image_height, image_width);
//
// 	// Write a tiff image.
// 	writeTiff_bw16(output_file_path_02, slice, image_height, image_width);
	
	// Free the array
	free(B);
	free(x);
	free(y);
	free(dp);
	free(I);
	
	// GTFO
	return(0);
}

























