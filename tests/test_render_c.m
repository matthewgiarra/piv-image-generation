% Image size
image_width = 1024;
image_height = 512;

num_particles = 1E5;

X = image_width * rand(num_particles, 1);
Y = image_height * rand(num_particles, 1);

dp = 3 * ones(size(X));
Io = 1 * ones(size(X));

% for k = 1 : num_particles
%     fprintf('X = %0.2f, Y = %0.2f\n', X(k), Y(k));
% end

tic;
B = generateParticleImage(image_height, image_width, X, Y, dp, Io);
toc;
% imagesc(B); axis image;
