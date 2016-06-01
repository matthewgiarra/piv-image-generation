#include <tiffio.h>
#include <string.h>
#include <stdint.h>

// Code for testing the C TIFF library, libtiff
void writeTiff_bw16(char *output_file_path, uint16_t *image_data, int image_height, int image_width){
		
	// Set the number of channels in the image to one, since we are assuming these images come out grayscale.
	const int samples_per_pixel = 1;
	
	// Counter
	int row;
						
	// Open the file
	TIFF *output_image = TIFFOpen(output_file_path, "w");
	
	/* Set some fields for the TIFF image */
	// Image width
	TIFFSetField(output_image, TIFFTAG_IMAGEWIDTH, image_width);
	
	// Image height
	TIFFSetField(output_image, TIFFTAG_IMAGELENGTH, image_height);
	
	// Samples per pixel
	TIFFSetField(output_image, TIFFTAG_SAMPLESPERPIXEL, samples_per_pixel);
	
	// Bits per sample. Set to 16 bits for these images.
	TIFFSetField(output_image, TIFFTAG_BITSPERSAMPLE, 16);
	
	// Set the origin of the image to the top-left corner.
	TIFFSetField(output_image, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
	
	// No compression
	TIFFSetField(output_image, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	
	// Some other important fields that I don't understand yet.
	TIFFSetField(output_image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(output_image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	
	// Set the length in memory of one row of pixels in the image.
	tsize_t line_bytes = samples_per_pixel * image_width;
		
	// Buffer used to store the row of pixel information for writing to a file
	uint16_t *buf = NULL;
	
	// Allocate memory to store the pixels of the current row
	if (TIFFScanlineSize(output_image) == line_bytes)
		// buf = (unsigned char*)_TIFFmalloc(line_bytes);
		buf = (uint16_t*)_TIFFmalloc(line_bytes);
	else
		buf = (uint16_t*)_TIFFmalloc(TIFFScanlineSize(output_image));
	
	// Set the strip size of the file to be the size of one row of pixels
	TIFFSetField(output_image, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(output_image, image_width * samples_per_pixel));
	
	// Write the image to the file one strip at a time
	for(row = 0; row < image_height; row++){
		memcpy(buf, &image_data[(row) * line_bytes], sizeof(uint16_t) * line_bytes);
		if(TIFFWriteScanline(output_image, buf, row, 0) < 0)
			break;
	}
	
	// Destroy the buffer
	_TIFFfree(buf);
	
	// Close the output image
	TIFFClose(output_image);
		
}













