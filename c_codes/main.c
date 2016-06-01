#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include "piv-image-gen.h"
#include "image_writing.h"

// Main function
int main(int argc, char *argv[]){

	if(argc < 6){
		printf("Usage: ./main <num_rows> <num_cols> <particles_per_pix> <num_threads> <output_path>\n");
		return(-1);
	}
	
	// Command line arguments
	const int num_rows = atoi(argv[1]);
	const int num_cols = atoi(argv[2]);
	float particles_per_pixel = atof(argv[3]);
	int max_num_threads = atoi(argv[4]);
	char *output_path = argv[5];
	
	// Number of particles
	int num_particles;
		
	// Allocate the particle coordinates
	float *X, *Y, *particle_diameters, *particle_max_intensities, *output_image;
	
	// Number of pixels in the image
	int num_pixels = num_rows * num_cols;
		
	// Counter	
	int p;
	
	// Number of particles
	num_particles = (int)round((float)num_pixels * (float)particles_per_pixel);
	
	// Number of particles
	printf("Number of particles: %d\n", num_particles);
	
	// Upper and lower bounds
	float lower_bound_x = 0;
	float upper_bound_x = num_cols;
	float lower_bound_y = 0;
	float upper_bound_y = num_rows;
	
	float lower_bound_dp = 1.8;
	float upper_bound_dp = 4.40;
	
	// Allocate the array
	X = malloc(num_particles * sizeof(float));
	Y = malloc(num_particles * sizeof(float));
	particle_diameters = malloc(num_particles * sizeof(float));
	particle_max_intensities = malloc(num_particles * sizeof(float));

	// Allocate the output image
	output_image = calloc(num_pixels, sizeof(float));
	
	// Populate the arrays
	devrand(X, num_particles, lower_bound_x, upper_bound_x);
	devrand(Y, num_particles, lower_bound_y, upper_bound_y);
	devrand(particle_diameters, num_particles, lower_bound_dp, upper_bound_dp);
	
	// Particle diameters
	for(p = 0; p < num_particles; p++){
		particle_max_intensities[p] = lower_bound_dp / particle_diameters[p];
	}

	// Generate the image
	generateParticleImage_omp(output_image, num_rows, 
		num_cols, num_particles, 
		X, Y, particle_diameters, 
		particle_max_intensities, max_num_threads);
		
	// Write the image
	printf("Image save path: %s\n", output_path);
			
	
	// Free the array
	free(X);
	free(Y);
	free(output_image);
	free(particle_diameters);
	free(particle_max_intensities);
	return(0);
}





