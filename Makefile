# Makefile for Gudbrand Tandberg's NBODY Project
# FYS3150 - Computational Physics fall 2014

CC=g++
C_FLAGS = -Wall -O3


SRCDIR = source
OBJDIR = objects

#SOURCES := $(wildcard $(SRCDIR)/*.cpp)
#INCLUDES := $(wildcard $(SRCDIR)/*.h)
#OBJS := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

MAIN_OBS = objects/main.o objects/Body.o objects/NBodySolver.o
GL_OBS = imageloader.o openGL_solarsystem3.o


#link armadillo:
ARMA_FLAGS = -larmadillo -framework Accelerate

#link GLUT and OpenGL:
GL_FLAGS = -framework GLUT -framework OpenGL


## rules for executables:

main: objects/main.o objects/Body.o objects/NBodySolver.o
	$(CC) $(C_FLAGS) -o $@ $^ $(ARMA_FLAGS)

openGL_3: objects/imageloader.o objects/openGL_solarsystem3.o
	$(CC) $^ -o $@ $(GL_FLAGS)

## rules for objects:
$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CC) -c $< -o $@

clean:
	rm -f main openGL_3 objects/*
