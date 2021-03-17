# LorenzSystem-RungeKutta
A Fortran 90 project to analyze the Lorenz System using the 4th order Runge-Kutta Method. It also includes plots of the phase space, the parameter space and the Lyapunov coefficient made with Wolfram Mathematica.

The lorenz.f90 program solves the Lorenz System using the 4th Order Runge-Kutta method and saves the data for plotting the Lorenz Map, the time series of the trajectories x, y, and z and the 2D and 3D coordinates of the phase space.

The lorenz2.f90 program solves the Lorenz System using the 4th Order Runge-Kutta method for a system with a small perturbation, then, using a least squares method obtains the Lyapunov coefficient and the data for the log-log plot of the displacement between a solution and it's perturbation. It also saves the data for plotting the time series of the trajectories x, y, and z and the 2D and 3D coordinates of the phase space.

Finally, the lorenz3.f90 solves the Lorenz System using the 4th Order Runge-Kutta method for every value of rho from 0 to 250 in steps of 0.25. Doing this we save the values of the maxima of the trajectories and obtain the bifurcation data. 

The Mathematica notebook has the plots of all the programs named above and the TeX folder includes a presentation with a brief description of what was done on the project.
