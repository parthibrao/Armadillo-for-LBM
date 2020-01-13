# Armadillo-for-LBM
Demonstration of the use of the Armadillo C++ linear algebra library for prototyping purposes in the lattice Boltzmann method. The repo contains code for single and two phase flows using BGK and MRT collision models. 
The directory consists of the following files:
1. ArmaFeatures.cpp: A quick rundown on the essential features of the Armadillo library in a LBM context
2. spinodalDecomp.cpp: Simulation of the 2D spinodal decompostion process using the Shan-Chen multicomponent LB model. The code illustrates how Armadillo's datastructures and associated functionalities that can be useful for LBM purposes
3. 2DPoiseulleMRT.cpp: Simulation of 2D poiseuille flow using the MRT collision model and explicit (Guo) forcing scheme
4. 3D_LDC.cpp: 3D lid driven cavity using the BGK model
5. twoSpheres.dat: a sample input file containing indicator values of a particular node (lattice site), i.e., if a particular node is a solid (=1), fluid (=0) or a boundary (=2) node.
6. plotSpinodalDecomp.py: A Python (matplotlib) file used for plotting results obtained from spinodalDecom.cpp
7. staticDropletBGK.cpp: Simulate droplet equilibrium using ShanChen twophase model + BGK collision + Guo (explicit) forcing
8. staticDropletMRT.cpp: Simulate droplet equilibrium using ShanChen twophase model + MRT collision + Guo (explicit) forcing


For more details on Armadillo API, refer to its documentation at http://arma.sourceforge.net/docs.html and the following paper:
Conrad Sanderson and Ryan Curtin. Armadillo: a template‐based C++ library for linear algebra. Journal of Open Source Software, Vol. 1, pp. 26, 2016.
