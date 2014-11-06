/*
 * openGL_solarsystem.cpp
 *
 * Part of Gravitational N-Body Simulations by Gudbrand Tandberg
 * FYS3150 Fall 2014
 *
 * Animates the sun-earth-moon system using OpenGL.
 *
 * USAGE:
 *    "r" key to toggle (off and on) running the animation
 *    "s" key to single-step animation
 *	  ESCAPE to exit.
 *	  "w", "s" to move camera forward/back
 *	  "left", "right" to rotate camera towards left/right
 */

#include "imageloader.h"
#include "openGL_solarsystem.h"
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<GLUT/glut.h>
#include<vector>
#include<cmath>

using namespace std;

// variables that control the animation state.

static GLenum spinMode = GL_TRUE;
static GLenum singleStep = GL_FALSE;

// number of planets to draw
int N = 100;
ifstream infile("/Users/gudbrand/Documents/NBODY/output/100_body_trajectories_100_0.1_0.dat");

// temporary coordinates
float x = 0;
float y = 0;
float z = 0;
vector<vector<float> > coordinates = vector<vector<float> >(N);

// variables for controlling the camera
static GLfloat camPos[3]={0,0,40};
static GLfloat lookAt[3]={0,0,0};

// to draw a texture on
GLUquadric *quad;

/*
 * glutKeyboardFunc is called below to set this function to handle
 * all normal key presses.
 */

 static void KeyPressFunc( unsigned char Key, int x, int y )
{
	switch (Key) {
		case 'R':
		case 'r':
			if ( singleStep ) {			// If ending single step mode
				singleStep = GL_FALSE;
				spinMode = GL_TRUE;		// Restart animation
			}
			else {
				spinMode = !spinMode;	// Toggle animation on and off.
			}
			break;
		case 'q':
		case 'Q':
			singleStep = GL_TRUE;
			spinMode = GL_TRUE;
			break;
			
		case 'w':
			camPos[0] += 0.1*(lookAt[0] - camPos[0]);
			camPos[2] += 0.1*(lookAt[2] - camPos[2]);
			break;
		case 'a':
			
			break;
		case 's':
			camPos[0] -= 0.1*(lookAt[0] - camPos[0]);
			camPos[2] -= 0.1*(lookAt[2] - camPos[2]);
			break;
		case 'd':
			
			break;
			
		case 27:	// Escape key
			exit(1);
	}
}

/* 
 * glutSpecialFunc is called below to set this function to handle
 * all special key presses.
 */

static void SpecialKeyFunc(int Key, int x, int y)
{
	switch (Key) {
			
		case GLUT_KEY_RIGHT:
			lookAt[0] = camPos[0] + cos(0.1)*(lookAt[0] - camPos[0]) - sin(0.1)*(lookAt[2] - camPos[2]);
			lookAt[2] = camPos[2] + sin(0.1)*(lookAt[0] - camPos[0]) + cos(0.1)*(lookAt[2] - camPos[2]);
			
			break;
		case GLUT_KEY_LEFT:
			lookAt[0] = camPos[0] + cos(-0.1)*(lookAt[0] - camPos[0]) - sin(-0.1)*(lookAt[2] - camPos[2]);
			lookAt[2] = camPos[2] + sin(-0.1)*(lookAt[0] - camPos[0]) + cos(-0.1)*(lookAt[2] - camPos[2]);
			break;
	}
}


/*
 * drawSphere() draws a textured sphere of given radius, where int planet 
 * tells which texture to use.
 */

static void drawSphere(int planet)
{
	float radius = 0.2;
	gluSphere(quad, radius, 20, 20);
}


/*
 * Animate() handles the animation and the redrawing of the
 * graphics window contents. Gets new coordinates from file
 * and draws spheres at these coordinates.
 */

static void Animate(void)
{
	// Clear the redering window
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (spinMode) {
		// Update the animation state (get new coordinates)
		
		for (int i=0; i<N; i++) {
			if (infile >> coordinates[i][0] >> coordinates[i][1] >> coordinates[i][2]) {
				
			}
			else{
				exit(1);
			}
		}
	}

	// Clear the current matrix (Modelview)
    glLoadIdentity();
	
	glMatrixMode(GL_PROJECTION);
	gluPerspective(60.0, (GLfloat) 1080./720, 0.1, 200);
	
	//set the camera position
	gluLookAt(camPos[0],camPos[1],camPos[2],
			  lookAt[0], lookAt[1], lookAt[2],
			  0, 1, 0);

	// Draw the bodies
	for (int i=0; i<N; i++){
		
		x = coordinates[i][0];
		y = coordinates[i][1];
		z = coordinates[i][2];
		
		glTranslatef(y, z, x);
		glRotatef(90, 1, 0, 0);
		drawSphere(i);
		glRotatef(-90, 1, 0, 0);
		glTranslatef(-y, -z, -x);
	}

	// Flush the pipeline, and swap the buffers
    glFlush();
    glutSwapBuffers();

	if (singleStep) {
		spinMode = GL_FALSE;
	}

	glutPostRedisplay();		// Request a re-draw for animation purposes

}


// ResizeWindow is called when the window is resized
static void ResizeWindow(int w, int h)
{
    float aspectRatio;
	h = (h == 0) ? 1 : h;
	w = (w == 0) ? 1 : w;
	glViewport( 0, 0, w, h );	// View port uses whole window
	aspectRatio = (float)w/(float)h;

	// Set up the projection view matrix (not very well!)
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective( 60.0, aspectRatio, 1.0, 30.0 );

	// Select the Modelview matrix
    glMatrixMode(GL_MODELVIEW);
}

// Load texture
GLuint loadTexture(Image* image) {
	GLuint textureId;
	glGenTextures(1, &textureId); //Make room for our texture
	glBindTexture(GL_TEXTURE_2D, textureId); //Tell OpenGL which texture to edit
	//Map the image to the texture
	glTexImage2D(GL_TEXTURE_2D,                //Always GL_TEXTURE_2D
				 0,                            //0 for now
				 GL_RGB,                       //Format OpenGL uses for image
				 image->width, image->height,  //Width and height
				 0,                            //The border of the image
				 GL_RGB, //GL_RGB, because pixels are stored in RGB format
				 GL_UNSIGNED_BYTE, //GL_UNSIGNED_BYTE, because pixels are stored
				 //as unsigned numbers
				 image->pixels);               //The actual pixel data
	return textureId; //Returns the id of the texture
}

// Initialize OpenGL's rendering modes and load planet textures
void OpenGLInit(void)
{
	
	glShadeModel( GL_FLAT );
	glClearColor( 0.0, 0.0, 0.0, 0.0 );
	glClearDepth( 1.0 );
	glEnable( GL_DEPTH_TEST );
	glEnable(GL_NORMALIZE);
	glEnable(GL_COLOR_MATERIAL);
	
	// set lighting
	GLfloat light_ambient[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat light_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_position0[] = { 0.0, 30, 0.0, 0.0 };

	
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position0);
	
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	
	//initialize coordinates
	for (int i=0; i<N; i++){
		coordinates[i] = vector<float>(3);
		for (int j=0; j<3; j++) {
			coordinates[i][j] = 0;
		}
	}
	
	quad = gluNewQuadric();
	
}

/*  Main program
 *  Set up OpenGL, hook up callbacks, and start the main loop
 */

 int main( int argc, char** argv )
{
	// Need to double buffer for animation
	glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );

	// Create and position the graphics window
    glutInitWindowPosition( 0, 0 );
    glutInitWindowSize( 1080, 720 );
    glutCreateWindow( "Solar System Demo" );

	// Initialize OpenGL.
    OpenGLInit();
	
	

	// Set up callback functions for key presses
	glutKeyboardFunc(KeyPressFunc);
	glutSpecialFunc(SpecialKeyFunc);

	// Set up the callback function for resizing windows
    glutReshapeFunc(ResizeWindow);
	
	// Callback for graphics image redrawing
    glutDisplayFunc(Animate);
	
	// Start the main loop.  glutMainLoop never returns.
	glutMainLoop();

    return(0);
}
