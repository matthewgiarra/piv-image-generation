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
void poiseuille(double *x, double *y, double *z, double *x_new, double yc, double zc, double h, double v_max, double num_points){
	
	// Allocate normalized channel position
	double r_squared;
	
	// Loop over all the points
	for(int k = 0; k < num_points; k++){
		
		// Coordinates centered about the origin.
		double y_cent = y[k] - yc;
		double z_cent = z[k] - zc;
		
		// Radial position
		r_squared = y_cent * y_cent + z_cent * z_cent;
		
		// Poiseuille equation
		x_new[k] = x[k] + v_max * (1 - r_squared / (h * h));
	}		
}