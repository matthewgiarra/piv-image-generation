x = X_PIX(inds);
y = Y_PIX(inds);
ints = min(5 * particle_max_intensities, 1);
dp = particle_image_diameters_pixels;
marker_areas = pi / 4 * dp.^2;

num_particles = length(x);

hold off

for k = 1 : num_particles
   scatter(x(k), y(k), marker_areas(k), 'ok', 'filled',...
       'markerfacealpha', ints(k), ...
       'markeredgecolor', 'none' );
%        
    
 hold on   
    
end

       
    


hold off
axis image;
ylim([1, 127]);
xlim([1, 127]);
