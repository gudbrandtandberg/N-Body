/*
 * imageloader.h
 *
 * Part of Gravitational N-body simulation project by Gudbrand Tandberg
 * FYS3150 fall 2014
 *
 * Third party program originally written by Bill Jacobs at videotutorialsrock.com
 * that loads a bitmap image into a Image object.
 */

#ifndef IMAGE_LOADER_H_INCLUDED
#define IMAGE_LOADER_H_INCLUDED

//Represents an image
class Image {
	public:
		Image(char* ps, int w, int h);
		~Image();

		/* An array of the form (R1, G1, B1, R2, G2, B2, ...) indicating the
		 * color of each pixel in image.  Color components range from 0 to 255.
		 * The array starts the bottom-left pixel, then moves right to the end
		 * of the row, then moves up to the next column, and so on.  This is the
		 * format in which OpenGL likes images.
		 */
		char* pixels;
		int width;
		int height;
};

//Reads a bitmap image from file.
Image* loadBMP(const char* filename);

#endif