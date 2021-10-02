# MMLBM
Multiblock and multigrid lattice Boltzmann methods in domain of complex geometries

This MATLAB library implements the examples of lattice Boltzmann methods in the thesis <br/>
[The Lattice Boltzmann Method for Fluid Dynamics: Theory and Applications](https://github.com/cpempire/MMLBM/tree/master/thesis)

<br/>

which includes 
- Poiseuille flow, 
- lid driven cavity flow, 
- Womersley flow, 
- Taylor-Couette flow, 
- blood flow in aneurysms, etc. 

<br/>
This library also includes the most common schemes for imposing boundary conditions on both straight and curved boundaries. Moreover, it features multiblock and multigrid implementation that enables the lattice Boltzmann methods to more practical applications, where multiscale physical problem in complex geometries.

<img src="images/cerebral_aneurysm_2.jpeg" style="width: 24%; max-width: 250px; margin-top: 1.5em;">
<img src="images/lbm-block.png" style="width: 24%; max-width: 250px; margin-top: 1.5em;">

<img src="images/lbm-velocity.png" style="width: 24%; max-width: 250px; margin-top: 1.5em;">
<img src="images/lbm-pressure.png" style="width: 24%; max-width: 250px; margin-top: 1.5em;">

![arm](images/lbm-plaque-velocity.png)

![arm](images/lbm-plaque-pressure.png)
