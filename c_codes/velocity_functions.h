void plane_poiseuille(double *x, double *y, double *x_new, double h, double v_max, double num_points){
	
	// Allocate normalized channel position
	double yh;
	
	// Loop over all the points
	for(int k = 0; k < num_points; k++){
		
		// Position in channel
		yh = y[k] / h;
		
		// Poiseuille equation
		x_new[k] = x[k] + v_max * (1 - yh * yh);
	}
		
}

// Poiseuille flow
void poiseuille(double *x, double *y, double *z, double *x_new, double h, double v_max, double num_points){
	
	// Allocate normalized channel position
	double r_squared;
	
	// Loop over all the points
	for(int k = 0; k < num_points; k++){
		
		// Radial position
		r_squared = y[k] * y[k] + z[k] * z[k];
		
		// Poiseuille equation
		x_new[k] = x[k] + v_max * (1 - r_squared);
	}		
}