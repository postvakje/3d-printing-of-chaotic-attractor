# 3d-printing-of-chaotic-attractor
Python port of the MATLAB file `runrucklidge.m` from "Modeling Dynamical Systems for 3D Printing”, Notices of the AMS, 2020 by Stephen K. Lucas, Evelyn Sander, and Laura Taalman and available at http://math.gmu.edu/~sander/EvelynSite/supplementary-materials-for.html
     
    
Some features added in the port: 
1. dynamical system descriptions are encapsulated in `dynamical_system` base class.
2. arclength extension to the system equation are in `dynamical_system` base class, i.e. in the definition of the dynamical system inherited from the base class only the original system equations need to be defined.
3. use scipy `minimize` module rather than the secant method to find the time when a certain length of the trajectory curve is reached.

The main file is `run_dynamical_system.py`. Put both files `run_dynamical_system.py` and `dynamical_systems.py` in the same directory. Run the file `run_dynamical_system.py` to generate an output file of data points. 4 dynamical systems are defined in `dynamical_systems.py`. The hyperchaotic Rossler system is a 4-dimensional dynamical system, but only the first 3 coordinates of the state vector are saved to the data file. A periodically driven nonautonomous chaotic system `nonauto_chaotic_system` is also included.

To switch dynamical systems, change the line `system=nonauto_chaotic_system()` in `run_dynamical_system.py` to another system.

![nonautonomous chaotic attractor](nonautonomous_chaotic_attractor.jpg)