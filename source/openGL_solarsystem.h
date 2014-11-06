/*
 * openGL_solarsystem.h
 *
 * Part of Gravitational N-Body Simulations by Gudbrand Tandberg
 * FYS3150 Fall 2014
 *
 * Animates the trajectories of the inner solar system as
 * calculated by NBodySolver with initial conditions as of
 * 1.10.2014
 *
 * Based on Solar.c by Samuel R. Buss taken from
 * http://www.math.ucsd.edu/~sbuss/MathCG/OpenGLsoft/Solar/Solar.c
 *
 *
 * USAGE:
 *    "r" key to toggle (off and on) running the animation
 *    "s" key to single-step animation
 *	  ESCAPE to exit.
 *	  "w", "s" to move camera forward/back
 *	  "left", "right" to rotate camera towards left/right
 */

#include "imageloader.h"
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<GLUT/glut.h>
#include<vector>
#include<cmath>

using namespace std;

/*
 *  Custom initialization of OpenGL. Sets up rendering modes and loads planet textures
 */

void OpenGLInit(void);

/*
 * Draws a textured sphere of given radius, at the surrent positions, 
 * where int planet tells which texture to use.
 */

static void drawSphere(int planet);

/*
 * Handles the animation and the redrawing of the
 * graphics window contents. Gets new coordinates from file
 * and draws spheres at these coordinates.
 */

static void Animate(void);

/*
 * Takes an Image and returns the id of this texture (a unsigned int)
 * written by Bill Jacobs at
 * http://www.videotutorialsrock.com/opengl_tutorial/reference.php#lighting
 */

GLuint loadTexture(Image* image);

/*
 * Handles resizing of the graphics window
 */

static void ResizeWindow(int w, int h);

/*
 * glutKeyboardFunc is called in main to set this function to handle
 * all normal key presses. Handles w, s, q and r keypresses.
 */

static void KeyPressFunc(unsigned char Key, int x, int y);

/*
 * glutSpecialFunc is called in main to set this function to handle
 * all special key presses. Handles left and right key presses.
 */

static void SpecialKeyFunc(int Key, int x, int y);

/*  
 * Main program
 * Set up OpenGL, hook up callbacks, and start the mainloop
 */

int main(int argc, char ** argv);


