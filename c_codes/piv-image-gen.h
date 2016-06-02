// Include the standard libraries.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Constants
#define PI 3.141592653589793

void calculate_particle_image_size(float *particle_image_diameters_pixels, float *particle_max_intensities, // Result arrays
float *Z_microns, float magnification, float particle_diameter_microns, float wavelength_microns, float NA, float pixel_size_microns, int num_particles);



void generateParticleImage(float *output_image, int num_rows, int num_cols, 
						int num_particles, float *X, float *Y, float *particle_diameters,
						float *particle_max_intensities, float dp_max) /* Input variables */

{

	
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
		maxRenderedCol = (int)ceil(X[p]  + render_fraction * dp[p]);
		minRenderedRow = (int)floor(Y[p] - render_fraction * dp[p]);
		maxRenderedRow = (int)ceil(Y[p]  + render_fraction * dp[p]);
    
        // Flag whether or not to render the particle
        render_particle = (minRenderedCol <= (num_cols - 1)) & (maxRenderedCol >= 0) &
                          (minRenderedRow <= (num_rows - 1)) & (maxRenderedRow >= 0) &
					  													(dp[p] <= dp_max);
	
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
						output_image[ind] += Io[p] * dp[p] * dp[p] * PI / 32 *
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

// This function seeds the random number generator.
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

// Generate micro PIV image.
void generate_micro_piv_image( float *output_image, float *X_PIX, float *Y_PIX, float *Z_microns, int num_particles,
								int num_rows, int num_cols, float particle_diameter_microns, 
								float pixel_size_microns, float particle_volume_fraction, 
								float channel_depth_microns, float objective_magnification,
								float NA, float working_distance_microns, float focal_length_microns,
								float wavelength_microns){
									
	// Particle image diameters
	float * particle_image_diameters_pixels;
	float * particle_max_intensities;
	
	// Allocate image diameters
	particle_image_diameters_pixels = malloc(num_particles * sizeof(float));
	
	// Allocate max intensities
	particle_max_intensities = malloc(num_particles * sizeof(float));
	
	
	// Limit of the near-field in microns
	// This is calculated from the description 
	// in Eckstein & Vlachos 2009 "Digital particle image velocimetry (DPIV)
	// robust phase correlation"
	// by setting the size of a particle (in microns) at 
	// a distance z from the focal plane (taken from equation 21) 
	// equal to twice the average particle spacing (concentration in particles /
	// / microns^3) ^(-1/3)//

	// Particle concentration in particles per cubic micron
	float particle_concentration = 6.0 * particle_volume_fraction / 
	    (PI * pow(particle_diameter_microns, 3));
	
	// Largest diameter of a particle that is considered "in focus," in pixels.
	float dp_max_pixels = pow(particle_concentration, -1.0/3.0) * objective_magnification / pixel_size_microns;
	
	// Calculate the particle image diameters and intensities								
	calculate_particle_image_size(particle_image_diameters_pixels, particle_max_intensities, // Result arrays
			Z_microns, objective_magnification, particle_diameter_microns, 
			wavelength_microns, NA, pixel_size_microns, num_particles);
			
	// Generate the image
	generateParticleImage(output_image, num_rows, num_cols, 
							num_particles, X_PIX, Y_PIX, particle_image_diameters_pixels,
							particle_max_intensities, dp_max_pixels); /* Input variables */
			
	free(particle_image_diameters_pixels);	
	free(particle_max_intensities);						
}

// Calculate particle image diameter in diffraction limited optics.
// Microscope model taken from Olsen // Adrian 2000
void calculate_particle_image_size(float *particle_image_diameters_pixels, float *particle_max_intensities, // Result arrays
float *Z_microns, float magnification, float particle_diameter_microns, float wavelength_microns, float NA, float pixel_size_microns, int num_particles){
	
	// Redefine some variables for readability
	float M2 = magnification * magnification;
	float dp2 = particle_diameter_microns * particle_diameter_microns;
	float NA2 = NA * NA;
	float L = wavelength_microns;
	float particle_image_diameter_microns;
	
	// Counter
	int n;
							  	
	// Point response diameter
	float ds = 2.44 * (magnification + 1) * L / (2 * NA);
	float ds2 = ds * ds;
  									
	// In focus diameter
	float de_focused = sqrt(M2 * dp2 + ds2);

	// Image diameters in microns
	// and max intensities.
	for(n = 0; n < num_particles; n++){
		particle_image_diameter_microns = sqrt( M2 * dp2 + ds2 + 4 * NA2 * M2 * Z_microns[n] * Z_microns[n]);
		particle_image_diameters_pixels[n] = particle_image_diameter_microns * magnification / pixel_size_microns;
		particle_max_intensities[n] = de_focused * de_focused / (particle_image_diameter_microns * particle_image_diameter_microns);	
	}								
}



// Convert (row, column) subscript to linear index
int sub2ind(int row, int col, int num_cols){
	return col + num_cols * row;
	
}