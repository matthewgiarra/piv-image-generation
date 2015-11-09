
JobList = MonteCarloImageGenerationJobFile_ND;

case_name = JobList(1).CaseName;

[I1, I2] = generateMonteCarloImageSet_ND(JobList);

num_images = size(I2, 3);

image_dir = ['~/Desktop/' case_name];
if ~exist(image_dir, 'dir');
	mkdir(image_dir)
end



close all;
for k = 1 : num_images
	fprintf('Writing image %d of %d...\n', k, num_images)
	%
	im = (I2(:, :, k));
	
	image_name = [case_name '_' num2str(k, '%05d') '.tif'];

	image_path = fullfile(image_dir, image_name);

	imwrite(im, image_path, 'compression', 'none');

	imagesc(im);
	axis image;
% 	caxis([0, intmax('uint16')]);
	pause(0.01);
	
end