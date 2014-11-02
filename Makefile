# Standard makefile-mal laget av Gudbrand for C++ programmer i FYS3150

CC=g++
C_FLAGS = -Wall -O3 -g

PROG = main
OBJS = Body.o NBodySolver.o main.o

#hvis armadillo.h skal brukes:
LIB_FLAGS = -larmadillo -framework Accelerate


#regel for executables:
$(PROG): $(OBJS)
	$(CC) $(C_FLAGS) -o $@ $(OBJS) $(LIB_FLAGS)


openGL_solarsystem3: openGL_solarsystem3.cpp
	$(CC) $^ -o $@ -framework GLUT -framework OpenGL

#regel for objekter
Body.o: Body.cpp
	$(CC) $(C_FLAGS) -c Body.cpp

NBodySolver.o: NBodySolver.cpp
	$(CC) -c $^

main.o: main.cpp
	$(CC) -c $^

openGL_solarsystem3.o: openGL_solarsystem3.cpp
	$(CC) -c $^

clean:
	rm -f ${PROG} *.o
