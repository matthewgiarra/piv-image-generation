#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <string.h>
#include "piv-image-gen.h"
#include "image_writing.h"

// Main function
int main(int argc, char *argv[]){

	if(argc < 7){
		printf("Usage: ./main <num_rows> <num_cols>" 
			" <particles_per_pix>"
				" <num_images> <num_threads> <output_base>\n");
		return(-1);
	}
	
	// Command line arguments
	const int num_rows = atoi(argv[1]);
	const int num_cols = atoi(argv[2]);
	float particles_per_pixel = atof(argv[3]);
	int num_images = atoi(argv[4]);
	int max_num_threads = atoi(argv[5]);
	char *output_file_base = argv[6];
	
	// Length of the input string
	int len = strlen(argv[6]);
		
	// Number of particles
	int num_particles;
		
	// Allocate the particle coordinates
	float *X, *Y;
	float *particle_diameters, *particle_max_intensities;
	
	// Allocate the output image
	float *output_image;
	
	// Allocate the uint16 output image
	uint16_t *output_image_uint16;
	
	// Declare the output file path
	char *output_file_path;
	
	// Number of pixels in the image
	int num_pixels = num_rows * num_cols;
		
	// Counter	
	int p, n;
	
	// Number of particles
	num_particles = (int)round((float)num_pixels * (float)particles_per_pixel);
	
	// Number of particles
	printf("Number of particles: %d\n", num_particles);
	
	// Upper and lower bounds
	float lower_bound_x = 0;
	float upper_bound_x = num_cols;
	float lower_bound_y = 0;
	float upper_bound_y = num_rows;
	
	float lower_bound_dp = 2.8;
	float upper_bound_dp = 8.8;
	
	// Max intensity of any particle.
	float upper_bound_IO = pow(2, 15) - 1;
	float lower_bound_IO = 0.1 * upper_bound_IO;
	
	// Print number of images
	printf("Number of images: %d\n", num_images);
	
	// Seed the random number generator
	seed_random_number_generator();
	
	#pragma omp parallel private(X, Y, particle_max_intensities, particle_diameters, output_image, output_image_uint16, output_file_path) num_threads(max_num_threads)
	{
	
		// Allocate the array
		X = malloc(num_particles * sizeof(float));
		Y = malloc(num_particles * sizeof(float));
		particle_diameters = malloc(num_particles * sizeof(float));
		particle_max_intensities = malloc(num_particles * sizeof(float));

		// Allocate the output image
		output_image = calloc(num_pixels, sizeof(float));
		output_image_uint16 = malloc(num_pixels * sizeof(uint16_t));
		
		// Allocate output file path
		// Output file name stuff
		output_file_path = malloc((len + 20) * sizeof(char));
	
		#pragma omp for
		// Loop over images
		for(n = 0; n < num_images; n++){
	
			sprintf(output_file_path, "%s%05d.tiff", output_file_base, n);
			// printf("%s\n", output_file_path);
	
			// Populate the arrays
			devrand(X, num_particles, lower_bound_x, upper_bound_x);
			devrand(Y, num_particles, lower_bound_y, upper_bound_y);
			devrand(particle_max_intensities, num_particles, lower_bound_IO, upper_bound_IO);
			devrand(particle_diameters, num_particles, lower_bound_dp, upper_bound_dp);

			// // Particle diameters
			// for(p = 0; p < num_particles; p++){
			// 	particle_max_intensities[p] = IO_max * lower_bound_dp / particle_diameters[p];
			// }

			// Generate the image
			generateParticleImage(output_image, num_rows,
				num_cols, num_particles,
				X, Y, particle_diameters,
				particle_max_intensities);

			// Convert to uint16
			for(p = 0; p < num_pixels; p++){
				output_image_uint16[p] = (uint16_t)output_image[p];
			}

			// Write the image
			// printf("Image save path: %s\n", output_file_path);
			writeTiff_bw16(output_file_path, output_image_uint16 , num_rows, num_cols);

		}
	
		// Free the array
		free(X);
		free(Y);
		free(output_image);
		free(output_image_uint16);
		free(particle_diameters);
		free(particle_max_intensities);
		free(output_file_path);	
		
	}
	


	return(0);
}





