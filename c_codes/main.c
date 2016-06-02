#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <string.h>
#include "piv-image-gen.h"
#include "image_writing.h"

// Main function
int main(int argc, char *argv[]){

	if(argc < 6){
		printf("Usage: ./main <num_rows> <num_cols>" 
				" <num_images> <num_threads> <output_base>\n");
		return(-1);
	}
	
	// Command line arguments
	const int num_rows = atoi(argv[1]);
	const int num_cols = atoi(argv[2]);
	int num_images = atoi(argv[3]);
	int max_num_threads = atoi(argv[4]);
	char *output_file_base = argv[5];
	
	// Length of the input string
	int len = strlen(argv[5]);
		
	// Number of particles
	int num_particles;
		
	// Allocate the particle coordinates
	float *X, *Y, *Z;
	
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
	
	// Some constants
	float particle_diameter_microns = 0.1;
	float objective_magnification   = 60;
	float channel_depth_microns     = 100;
	float wavelength_microns        = 0.532;
	float pixel_size_microns 		= 10;
	float particle_volume_fraction  = 5E-5;
	float NA = 1.4;
	float working_distance_microns  = 0.3E3;
	float focal_length_microns      = 3.3333E3;
	float int_max = 0.95 * pow(2, 16) - 1;
	
	// Size of the domain
	float domain_width_microns  = (float)num_cols * pixel_size_microns / objective_magnification;
	float domain_height_microns = (float)num_rows * pixel_size_microns / objective_magnification;
	float domain_volume = domain_width_microns * domain_height_microns * channel_depth_microns;
	
	// Particle concentration in particles per cubic micron
	float particle_concentration = 6 * particle_volume_fraction / 
	    (PI * pow(particle_diameter_microns, 3));
	
	// Number of particles
	num_particles = (int)round(particle_concentration * domain_volume);
	
	// Number of particles
	printf("Number of particles: %d\n", num_particles);
	
	// Upper and lower bounds
	float lower_bound_x = 0;
	float upper_bound_x = num_cols;
	float lower_bound_y = 0;
	float upper_bound_y = num_rows;
	float lower_bound_z = -1 * channel_depth_microns / 2.0;
	float upper_bound_z =  1 * channel_depth_microns / 2.0;

	// Print number of images
	printf("Number of images: %d\n", num_images);
	
	// Seed the random number generator
	seed_random_number_generator();
	
	// #pragma omp parallel private(X, Y, Z, output_image, output_image_uint16, output_file_path) num_threads(max_num_threads)
	// {
	
		// Allocate the array
		X = malloc(num_particles * sizeof(float));
		Y = malloc(num_particles * sizeof(float));
		Z = malloc(num_particles * sizeof(float));

		// Allocate the output image
		output_image = calloc(num_pixels, sizeof(float));
		output_image_uint16 = malloc(num_pixels * sizeof(uint16_t));
		
		// Allocate output file path
		// Output file name stuff
		output_file_path = malloc((len + 20) * sizeof(char));
	
		// #pragma omp for
		// Loop over images
		for(n = 0; n < num_images; n++){
	
			sprintf(output_file_path, "%s%05d.tiff", output_file_base, n);
			// printf("%s\n", output_file_path);
	
			// Populate the arrays
			devrand(X, num_particles, lower_bound_x, upper_bound_x);
			devrand(Y, num_particles, lower_bound_y, upper_bound_y);
			devrand(Z, num_particles, lower_bound_z, upper_bound_z);
			
			// Generate the micro PIV image.
			generate_micro_piv_image(output_image, X, Y, Z, num_particles,
											num_rows, num_cols, particle_diameter_microns, 
											pixel_size_microns, particle_volume_fraction, 
											channel_depth_microns, objective_magnification,
											NA, working_distance_microns, focal_length_microns,
											wavelength_microns);
											
			// Pixel gain
			for(p = 0; p < num_pixels; p++){
				output_image[p] = int_max * output_image[p];
			}
			
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
		free(output_file_path);	
	// }
	
	return(0);
}





