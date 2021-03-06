# Makefile for Gudbrand Tandberg's NBODY Project
# FYS3150 - Computational Physics fall 2014

PCC = MPIC++
CC = /usr/bin/clang
C_FLAGS = -O3 -Wall

SRCDIR = source
OBJDIR = objects

#-lstdc++


MAIN_OBS = objects/main.o objects/Body.o objects/NBodySolver.o
SOLAR_OBS = objects/imageloader.o objects/openGL_solarsystem.o
CLUSTER_OBS = objects/imageloader.o objects/openGL_cluster.o

#link armadillo:
#-framework Accelerate
ARMA_FLAGS = -larmadillo

#link GLUT and OpenGL:
GL_FLAGS = -framework GLUT -framework OpenGL

OMP_FLAGS = -fopenmp

## rules for executables:

main: $(MAIN_OBS)
	$(CC) -v $(C_FLAGS) -o $@ $^ $(ARMA_FLAGS) 

solar_animation: $(SOLAR_OBS)
	$(CC) $(C_FLAGS) $^ -o $@ $(GL_FLAGS)

cluster_animation: $(CLUSTER_OBS)
	$(CC) $(C_FLAGS) $^ -o $@ $(GL_FLAGS)

## rules for objects:

$(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CC) $(C_FLAGS) -c $< -o $@ 

clean:
	rm -f main solar_animation cluster_animation objects/*
