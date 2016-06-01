// Include the standard libraries.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void generateParticleImage(float *output_image, int num_rows, 
	int num_cols, int num_particles, float *X, float *Y, float *particle_diameters, float *particle_max_intensities) /* Input variables */

{
	// Constants
    #define pi 3.141592653589793
	
	// Square root of 8
	#define sqrt8 2.82842712475
    
    // Machine precision
    #define eps 2.2204E-16
	
	// Render fraction
	#define render_fraction 0.75
	
	// Particle diameters
	#define dp particle_diameters
	
	// Max intensities
	#define Io particle_max_intensities
    		
	// Counters for loops. 
	// p is particle number;
	// r is row number;
	// c is column number;
	// ind is the pixel index 
	// within the output array (B)
	int p, r, c, ind; 
	
	// Radial coordinates from the
	// frame of reference of a particle
	double render_radius;
	
	// Integer rows and columns
	//to be rendered for each particle
	int minRenderedCol;
	int maxRenderedCol;
	int minRenderedRow;
	int maxRenderedRow;
	
    // Boolean render flag
    int render_particle, render_pixel;
		
	// Number of pixels
	int num_pixels = (num_rows) * (num_cols);
	
	// Set the output image to zeros
	memset(output_image, 0, num_pixels * sizeof(float));
		
	// Loop over all the particles		
	for(p = 0; p < num_particles; p++){		
	
		// Min and max rows and columns to render for that column.
		minRenderedCol = (int)floor(X[p] - render_fraction * dp[p]);
		maxRenderedCol = (int)ceil(X[p] + render_fraction * dp[p]);
		minRenderedRow = (int)floor(Y[p] - render_fraction * dp[p]);
		maxRenderedRow = (int)ceil(Y[p] + render_fraction * dp[p]);
    
        // Flag whether or not to render the particle
        render_particle = (minRenderedCol <= (num_cols - 1)) & (maxRenderedCol >= 0) &
                          (minRenderedRow <= (num_rows - 1)) & (maxRenderedRow >= 0);
	
		// Render the particle if the criteria are met.
        if(render_particle){
        
			// Increment the number of particles rendered
			// particles_rendered_private++;
		
			// Loop over all the particles.
			for(r = minRenderedRow; r <= maxRenderedRow; r++){
				for(c = minRenderedCol; c <= maxRenderedCol; c++){
				
					// Radius from particle center
					render_radius = sqrt(pow(c - X[p], 2) + pow(r - Y[p], 2));
				
					// Flag whether or not to render the pixel
					render_pixel = (c >= 0) & (c <= num_cols - 1) &
								   (r >= 0) & (r <= num_rows - 1) &
								   (render_radius <= render_fraction * dp[p]); 
				
					// Render the pixel if the criteria are met.
					if(render_pixel){
					
						// Index of the pixel to render					
						ind = r + (num_rows) * c;
				
						// Add the intensity to the image
						output_image[ind] += Io[p] * dp[p] * dp[p] * pi / 32 *
						                   (erf( sqrt8 * (c - X[p] - 0.5)
						                   / dp[p]) - erf(sqrt8 *
						                   (c - X[p] + 0.5) / dp[p])) *
						                   (erf( sqrt8 * (r - Y[p] - 0.5)
						                   / dp[p]) - erf(sqrt8 *
						                   (r - Y[p] + 0.5) / dp[p]));	
						}	
					}
				}	
	        }
		}
			
	 			
    return;	
}

// Random number generator.
void devrand(float *array, int array_length, float lower_bound, float upper_bound) {
	
	// Counter
	int k;
	
	// Generate the numbers between the upper and lower bounds.
	for (k = 0; k < array_length; k++) {
		array[k] = rand() / (float)RAND_MAX * (upper_bound - lower_bound) + lower_bound ;
	}
	
}

int seed_random_number_generator(){
	
	// Declare the seed
	int *seed;
	seed = malloc(sizeof(int));
	
	// Declare the urandom file pointer
	FILE *urandom;
	
	// Open the /urandom file for reading
	urandom = fopen("/dev/urandom", "r");
	
	// Inform the user if the file wasn't opened successfully.
	if (urandom == NULL) {
		printf("Cannot open /dev/urandom!\n");
		exit(-1);
	}
	
	// Put a random number from /urandom into the memory location of ths seed.
	fread(seed, sizeof(int), 1, urandom);
	
	// Seed the random number generator
	srand(*seed);
	
	// Close urandom
	fclose(urandom);
	
	// Free the seed
	free(seed);
	
	// Return 0 on success
	return(0);
}


// Convert (row, column) subscript to linear index
int sub2ind(int row, int col, int num_cols){
	return col + num_cols * row;
	
}