/*
 * openGL_solarsystem.cpp
 *
 * Part of Gravitational N-Body Simulations by Gudbrand Tandberg
 * FYS3150 Fall 2014
 *
 */

#include "openGL_solarsystem.h"

// variables that control the animation state.

static GLenum spinMode = GL_TRUE;
static GLenum singleStep = GL_FALSE;

// number of planets to draw
int N = 6;
ifstream infile("./output/6_body_trajectories_619_0.05_0_0.dat");

// temporary coordinates
float weekOfYear = 0;
float lengthOfDay[] = {25./7, 58./7, 243./7, 1./7, 1./7, 0.4/7};

float dt = 0.1;  //weeks
float x = 0;
float y = 0;
float z = 0;
vector<vector<float> > coordinates = vector<vector<float> >(N);

// variables for controlling the camera
static GLfloat camPos[3]={0,0.3,3.5};
static GLfloat lookAt[3]={0,0.3,0};

// texture for each body
GLuint _sunTexture;
GLuint _earthTexture;
GLuint _mercuryTexture;
GLuint _venusTexture;
GLuint _marsTexture;
GLuint _jupiterTexture;
GLuint _starsTexture;

// to draw a textured sphere
GLUquadric *quad;


 static void KeyPressFunc( unsigned char Key, int x, int y )
{
	switch (Key) {
		case 'R':
		case 'r':
			if (singleStep) {			// If ending single step mode
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
			camPos[0] += 0.2*(lookAt[0] - camPos[0]);
			camPos[2] += 0.2*(lookAt[2] - camPos[2]);
			break;
		case 'a':
			
			break;
		case 's':
			camPos[0] -= 0.2*(lookAt[0] - camPos[0]);
			camPos[2] -= 0.2*(lookAt[2] - camPos[2]);
			break;
		case 'd':
			
			break;
			
		case 27:	// Escape key
			exit(1);
	}
}

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

static void drawSphere(int planet)
{
	float radius;
	GLuint texture;
	switch (planet) {
		case 0:
			radius = 0.2;
			texture = _sunTexture;
			break;
		case 1:
			radius = 0.01;
			texture = _mercuryTexture;
			break;
		case 2:
			radius = 0.03;
			texture = _venusTexture;
			break;
		case 3:
			radius = 0.04;
			texture = _earthTexture;
			break;
		case 4:
			radius = 0.04;
			texture = _marsTexture;
			break;
		case 5:
			radius = 0.1;
			texture = _jupiterTexture;
			break;
		case 6:
			radius = 0.1;
			texture = _jupiterTexture;
			break;
		default:
			radius = 0.1;
			texture = _earthTexture;
			break;
	}
	
	//glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, texture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	gluQuadricTexture(quad, true);
	glRotatef(-weekOfYear/lengthOfDay[planet], 0, 0, 1);
	gluSphere(quad, radius, 20, 20);
	glRotatef(weekOfYear/lengthOfDay[planet], 0, 0, 1);
}

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
		weekOfYear += 1;
	}
	
	// Clear the current matrix
    glLoadIdentity();
	
	//set the camera position
	gluLookAt(camPos[0],camPos[1],camPos[2],
			  lookAt[0], lookAt[1], lookAt[2],
			  0, 1, 0);
	
	// Draw the planets
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

	glutPostRedisplay(); // Request a re-draw for animation purposes

}


static void ResizeWindow(int w, int h)
{
    float aspectRatio;
	h = (h == 0) ? 1 : h;
	w = (w == 0) ? 1 : w;
	glViewport(0, 0, w, h);	// View port uses whole window
	aspectRatio = (float)w/(float)h;

	// Set up the projection view matrix (not very well!)
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0, aspectRatio, 1.0, 30.0);

	// Select the Modelview matrix
    glMatrixMode(GL_MODELVIEW);
}


GLuint loadTexture(Image* image)
{
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
	return textureId;
}

void OpenGLInit(void)
{
	glShadeModel(GL_SMOOTH);
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClearDepth(1.0);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_NORMALIZE);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_TEXTURE_2D);

	// set lighting
	GLfloat light_ambient[] = {0.7, 0.7, 0.7, 1.0};
	GLfloat light_diffuse[] = {0.6, 0.6, 0.6, 1.0};
	GLfloat light_specular[] = {0.3, 0.3, 0.3, 0.3};
	GLfloat light_position0[] = {0.0, 0.0, 0.0, 0.0};

	
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

	//Load images & apply textures.
	
	Image* image = loadBMP("./images/sun.bmp");
	_sunTexture = loadTexture(image);
	
	image = loadBMP("./images/earth.bmp");
	_earthTexture = loadTexture(image);
	
	image = loadBMP("./images/venus.BMP3");
	_venusTexture = loadTexture(image);
	
	image = loadBMP("./images/mars.BMP3");
	_marsTexture = loadTexture(image);
	
	image = loadBMP("./images/mercury.BMP3");
	_mercuryTexture = loadTexture(image);
	
	image = loadBMP("./images/jupiter.BMP3");
	_jupiterTexture = loadTexture(image);
	
	delete image;
	
}

int main(int argc, char** argv)
{
	glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

	// Create and position the graphics window
    glutInitWindowPosition( 0, 0 );
    glutInitWindowSize( 1080, 720 );
    glutCreateWindow("Inner Solar System" );

	// custom initialization
	OpenGLInit();

	glutKeyboardFunc(KeyPressFunc);
	glutSpecialFunc(SpecialKeyFunc);
    glutReshapeFunc(ResizeWindow);
    glutDisplayFunc(Animate);
	glutMainLoop();

    return(0);
}
