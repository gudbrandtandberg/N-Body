/*
 * openGL_solarsystem.h
 *
 * Part of Gravitational N-Body Simulations by Gudbrand Tandberg
 * FYS3150 Fall 2014
 *
 * Animates the trajectories of the inner solar system as
 * calculated by NBodySolver with initial conditions as of
 * 1.10.2014
 */

void OpenGLInit(void);

static void drawSphere(int planet);
static void Animate(void);
static void ResizeWindow(int w, int h);
static void KeyPressFunc(unsigned char Key, int x, int y);
static void SpecialKeyFunc(int Key, int x, int y);

