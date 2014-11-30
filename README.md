# N-Body Gravitational simulations

Main page for Gudbrand Tandberg's N-Body simulations project for FYS3150 Computational Physics
at the University of Oslo fall 2014. 

## Project structure:

### Latex
Contains the project report and all resources associated with it.

### articles
A collection of articles on and around the subject matter. Some referenced from the paper. 

### images
Bitmap images of all planets, and a simple dot image used in matlab plotting.

### initial_conditions
Files containing initial conditions for various configurations. Uses format

mass x y z vx vy vz \n
...

### matlab
Matlab programs that analyze the data provided by main

### output
Folder containing all output from main. Mainly trajectories and energies. Trajectory format

x1 y1 z1 x2 y2 z2 ... xN yN zN \n
...

### source
All cpp header and source code files. Most notably main.cpp, NBodySolver.cpp and *_gl.cpp

