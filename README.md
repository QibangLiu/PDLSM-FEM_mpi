# PDLSM-FEM_mpi
PDLSM-FEM solver is an open-source parallel implementation of coupled peridynamics least squares minimization and finite element method (PDLSM-FEM) 
in 2D and 3D by MPI technique. This solver is written in a C++ environment on cross platforms (Windows, Linux). It includes implicit static, explicit
dynamic, and implicit dynamic solvers for structure analysis under displacement or traction loading. 

  - `Pre-process` - It requires a mesh data file with the specified format. Users may use commercial finite element codes such as ANSYS, ABAQUS to generate the mesh data file.  
  - `Solver` - PDLSM-FEM includes implicit and explicit solvers. For implicit solvers, it stores the large sparse matrix in a compressed sparse row (CSR) format and solvers the large linear systems of equations by the Parallel Direct Sparse Solver (PDSS) of [Intel Math Kernel Library](https://software.intel.com/content/www/us/en/develop/tools/oneapi/base-toolkit/download.html?operatingsystem=window&distributions=webdownload&options=offline).

  - `Post-process` - It stores the results in VTK format, which can be visualized directly by [Paraview](https://www.paraview.org/).

See the [manual](https://github.com/QibangLiu/PDLSM-FEM_mpi/blob/master/PDLSM_FEM_solver%20manual.pdf) file for further details.  

# Examples

| <img src="examples\DiagPlate\diagPlateUY.gif" width="400"> | <img src="examples\DiagPlate\diagPlateSY.gif" width="400"> | 
| :---: | :---: | 
| [Displacement UY of diagonal plate]| [Stress SY of diagonal plate]|

| <img src="examples\DCT\DCT_UY.gif" width="400"> | <img src="examples\DCT\DCT_SY.gif" width="400"> | 
| :---: | :---: | 
| [Displacement UY of disk-shaped compact tension ] | [Stress SY of disk-shaped compact tension]|


# Developers

Qibang Liu, qibangliu@ksu.edu

Dept. of Mechanical & Nuclear Engineering, Kansas State University.

# References
  Liu, Q., Xin, X. J., Ma, J., & Wang, Y. (2021). Simulating quasi-static crack propagation by coupled peridynamics least square minimization with finite element method. Engineering Fracture Mechanics, 252, 107862. https://doi.org/10.1016/j.engfracmech.2021.107862

  Liu, Q., Xin, X. J., & Ma, J. (2021). Coupled Peridynamics Least Square Minimization with Finite Element Method in 3D and Implicit Solutions by Message Passing Interface. Journal of Peridynamics and Nonlocal Modeling. https://doi.org/10.1007/s42102-021-00060-3



