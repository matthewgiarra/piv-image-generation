// whittaker_blackman7.c
#include <math.h>
#include <stdio.h>
#include "mex.h" /* Always include this */
void mexFunction(int nlhs, mxArray *plhs[], /* Output variables */
int nrhs, const mxArray *prhs[]) /* Input variables */
{
	
	// INPUTS
	
	// nlhs is the number of output (left-side) arguments, or the size of the plhs array.
	
	// plhs is the array of output arguments.
	
	// nrhs is the number of input (right hand side) arguments, or the size of the prhs array.
	
	// prhs is the array of input arguments.

    #define B_OUT plhs[0]

    #define HEIGHT_IN prhs[0]
	
    #define WIDTH_IN prhs[1]
	
	#define X_IN prhs[2]

    #define Y_IN prhs[3]
	
	#define DP_IN prhs[4]
	
	#define MAX_INT_IN prhs[5]
	
	// Pointers to arrays
	// B is the output image
	// dp is the vector of particle image diameters
	// max_int is the vector of particle max intensities
	// X is the vector of the particles' horizontal positions (non-integer)
	// Y is the vector of the particles' vertical positions (non-integer)
    double *B, *dp, *I, *X, *Y;
    const double eps = 1.0E-15;
    const double pi = 3.141592653589793;
	const double sqrt8 = sqrt(8);
	
	int num_particles;
	
	// Height and width of the image to be rendered
	int M = (int)*mxGetPr(HEIGHT_IN);
	int N = (int)*mxGetPr(WIDTH_IN);
	int p, r, c, ind; // counters
	
	// Number of rows in the position vector
	int particle_position_rows = mxGetM(X_IN);
	int particle_position_cols = mxGetN(X_IN);

	// Particle diameters
	dp = mxGetPr(DP_IN);
	
	// Particle intensities
	I = mxGetPr(MAX_INT_IN);
	
	/* Get the pointer to the data of X */
    X = mxGetPr(X_IN);
	
	/* Get the pointer to the data of Y */
    Y = mxGetPr(Y_IN); 
	
	/* Create the output matrix */
    B_OUT = mxCreateDoubleMatrix(M, N, 0); 
	
	/* Get the pointer to the data of B */
    B = mxGetPr(B_OUT); 
	
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
	
	// Count the number of particles.
	// Accept either row or column vectors for the positions.
	num_particles = fmax(particle_position_rows, particle_position_cols);
		
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







