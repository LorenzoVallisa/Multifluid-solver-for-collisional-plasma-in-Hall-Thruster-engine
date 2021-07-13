# Multi-fluid solver for collisional plasma simulation in a full electromagnetic interaction inside an Hall-Thurster engine

Solver for Euler set of equations for collisional ions and electrons in their full EM interaction. Electrons and ions are hence solved without viscous dissipation closure term,
each by means of momentum equations. Electric source term is then computed thourgh solution of Possion equation, whereas magnetic filed is constant and imposed.

Downlaod and compile but  NOTE THAT before that you need to install the LIS library (https://www.ssisc.org/lis/index.en.html) and compile it, the dynamic linking is already in the Makefile.

# PDE systems
################################################

Euler system has 4 equations:\
rho    --- density\
rho ux --- x-momentum\
rho uy --- y-momentum\
rho E  --- total energy

All coupled with resolution of electric field through the Possion equation for electrostatic potential

# GRID
################################################

The solution U has dimension U(N_EQ, N_cells), and includes the two ghost cells (GC).
Numbering of cells and interfaces is as follows:

   
         cell1     int1   cell2    int2   cell3   int3          intN     cellN
    ----(GC_1)------|------(C)------|------(C)-----|--- (...) ---|------(GC_N)----
                  x_min                                        x_max
    
N_cells     = total number of cells
N_cells - 2 = number of physical cells
N_int = N_cells - 1 = number of interfaces

First and last interfaces coincide with the limits of the domain, x_min, x_max

There is also the possibility to run simulation using a non-uniform grid (Chebyshev's grid points)

