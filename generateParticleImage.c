// Include the standard libraries.
#include <stdio.h>
#include <math.h>
#include "mex.h" /* Always include this */

void mexFunction(int nlhs, mxArray *plhs[], /* Output variables */
int nrhs, const mxArray *prhs[]) /* Input variables */

{
	#define B_OUT plhs[0] // Output image (array)
	#define M_IN  prhs[0] // Height of the output image (integer)
	#define N_IN  prhs[1] // Width of the output image (integer)
	#define X_IN  prhs[2] // Horizontal positions of particles (vector)
	#define	Y_IN  prhs[3] // Vertical positions of particles (vector)
	#define DP_IN prhs[4] // Particle diameters (vector)
	#define I_IN  prhs[5] // Particle max intensities
	
	// Constants
    #define pi 3.141592653589793
	
	// Square root of 8
	#define sqrt8 2.82842712475
    
    // Machine precision
    #define eps 2.2204E-16
    
    // Render fraction
    #define render_fraction 0.75
    
	// Declare variables
	// B is the output image;
	// X and Y are the horizontal and
	// vertical positions of the particles;
	// dp is the array of particle diameters;
	// I is the array of particle max intensities.
	// M and N are the numbers of rows
	// and columns in the image to be generated.
	double *B, *X, *Y, *dp, *I, *M, *N;
		
	// Number of particles in the image
	int num_particles;
	
	// Counters for loops. 
	// p is particle number;
	// r is row number;
	// c is column number;
	// ind is the pixel index 
	// within the output array (B)
	int p, r, c, ind; 
	
	// Fractional rows and columns
	//to be rendered for each particle
	double minRenderedCol_fract;
	double maxRenderedCol_fract;
	double minRenderedRow_fract;
	double maxRenderedRow_fract;
	double render_radius;
	
	// Integer rows and columns
	//to be rendered for each particle
	int minRenderedCol;
	int maxRenderedCol;
	int minRenderedRow;
	int maxRenderedRow;
	
    // Boolean render flag
    int render_particle, render_pixel;
    
	// Get the pointer to the vectors containing particle positions
	X = mxGetPr(X_IN);
	Y = mxGetPr(Y_IN);
	
	// Get pointers to particle diameters and intensities
	dp = mxGetPr(DP_IN);
	I  = mxGetPr(I_IN);
	
	// Get pointers to image sizes
	M = mxGetPr(M_IN);
	N = mxGetPr(N_IN);
	
	// Measure array sizes (number of particles)
	num_particles = mxGetM(X_IN);
	
	// Allocate the output array (real valued)
	B_OUT = mxCreateDoubleMatrix(*M, *N, mxREAL);
	
	// Get pointer to the data in B
	B = mxGetPr(B_OUT);
        
	// Print number of particles
	printf("Number of particles: %d\n", num_particles);
		
	// Loop over all the particles
	for(p = 0; p < num_particles; p++){		
		
		// Min and max rows and columns to render for that column.
		minRenderedCol = (int)floor(X[p] - render_fraction * dp[p]);
		maxRenderedCol = (int)ceil(X[p] + render_fraction * dp[p]);
		minRenderedRow = (int)floor(Y[p] - render_fraction * dp[p]);
		maxRenderedRow = (int)ceil(Y[p] + render_fraction * dp[p]);
        
        // Flag whether or not to render the particle
        render_particle = minRenderedCol <= (*N-1) & maxRenderedCol >= 0 &
                          minRenderedRow <= (*M-1) & maxRenderedRow >= 0;
		
		// Render the particle if the criteria are met.
        if(render_particle == true){
            
			// Inform the user
			printf("Rendering particle %d of %d\n", p, num_particles);
			
			// Loop over all the particles.
			for(r = minRenderedRow; r <= maxRenderedRow; r++){
				for(c = minRenderedCol; c <= maxRenderedCol; c++){
					
					// Radius from particle center
					render_radius = sqrt(pow(c - X[p], 2) + pow(r - Y[p], 2));
					
					// Flag whether or not to render the pixel
					render_pixel = c >= 0 & c <= *N - 1 &
								   r >= 0 & r <= *M - 1 &
								   render_radius <= render_fraction * dp[p]; 
					
					// Render the pixel if the criteria are met.
					if(render_pixel == true){
						
						// Index of the pixel to render					
						ind = r + (*M) * c;
					
						// Add the intensity to the image
						B[ind] += I[p] * dp[p] * dp[p] * pi / 32 *
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

