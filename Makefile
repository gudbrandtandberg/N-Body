# Standard makefile-mal laget av Gudbrand for C++ programmer i FYS3150

CC=g++
C_FLAGS = -Wall -O3 -g

PROG = main
OBJS = Body.o NBodySolver.o main.o NBody_functions.o

#hvis armadillo.h skal brukes:
LIB_FLAGS = -larmadillo -framework Accelerate 

# hvis engine.h skal brukes:
#LIB_FLAGS = -L/Applications/MATLAB_R2014a.app/bin/maci64 -lmx -leng $(LIB_FLAGS)
#C_FLAGS= -I/Applications/MATLAB_R2014a.app/extern/include $(C_FLAGS)

#regel for executables:
$(PROG): $(OBJS)
	$(CC) $(C_FLAGS) -o $@ $(OBJS) $(LIB_FLAGS)

#regel for objekter
Body.o: Body.cpp
	$(CC) $(C_FLAGS) -c Body.cpp

NBody_functions.o: NBody_functions.cpp
	$(CC) -c $^

NBodySolver.o: NBodySolver.cpp
	$(CC) -c $^

main.o: main.cpp
	$(CC) -c $^

clean:
	rm -f ${PROG} *.o
