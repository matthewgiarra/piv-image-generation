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
    #define eps 2.2204E-16;
	
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
	
	// Integer rows and columns
	//to be rendered for each particle
	int minRenderedCol;
	int maxRenderedCol;
	int minRenderedRow;
	int maxRenderedRow;
	
	// Cutoff intensity for rendering
	double cutoff_intensity;
		
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
	
	// Specify the cutoff intensity
	// for rendering Gaussian particles.
	// Pixels are rendered when
	// the intensity of the particle
	// at the pixel is greater than
	// or equal to this value.
	// This comment is wrong!
	// This is not what this variable does!
    cutoff_intensity = eps;
    
    // Temporary counter
    int pr = 0;
		
	// Loop over all the particles
	for(p = 0; p < num_particles; p++){		
		
		// Min and max rows and columns to render for that column.
		minRenderedCol_fract = X[p] - 0.75 * dp[p];
		maxRenderedCol_fract = X[p] + 0.75 * dp[p];
		minRenderedRow_fract = Y[p] - 0.75 * dp[p];
		maxRenderedRow_fract = Y[p] + 0.75 * dp[p];
		
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
		if(minRenderedRow >=0 & maxRenderedRow < (int)*M & minRenderedCol >=0 & maxRenderedCol < (int)*N & I[p] > cutoff_intensity){
			
            pr += 1;
            
            printf("Rendered %d particles\n", pr);
            
			// Loop over all the particles.
			for(r = minRenderedRow; r <= maxRenderedRow; r++){
				for(c = minRenderedCol; c <= maxRenderedCol; c++){
					
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
	
    return;	
}

