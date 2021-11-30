---
title: 'PDLSM-FEM: Solver of Coupled Peridynamics Least Squares Minimization with Finite Element Method'
tags:
  - Peridynamics
  - Finite element method
  - Least square minimization
  - MPI
  - C++
authors:
  - name: Qibang Liu
    orcid: 0000-0001-7935-7907
    affiliation: 1
  - name: X.J. Xin
    affiliation: 1
affiliations:
 - name: Department of Mechanical and Nuclear Engineering, Kansas State University, Manhattan, Kansas, USA
   index: 1
date: 22 August 2021
bibliography: paper.bib
---

# Summary

'PDLSM-FEM solver' is a parallel implementation of coupled peridynamics least squares minimization and finite element method (PDLSM-FEM) in 2D and 3D by MPI technique. This solver is written in a C++ environment on cross platforms. It includes implicit static, explicit dynamic, and implicit dynamic solvers for structure analysis under displacement or traction loading. PDLSM-FEM solver stores the large sparse matrix in a compressed sparse row (CSR) format rather than a full matrix format, and solves the large linear systems of equations by the Parallel Direct Sparse Solver (PDSS) of [Intel Math Kernel Library](https://software.intel.com/content/www/us/en/develop/tools/oneapi/base-toolkit/download.html?operatingsystem=window&distributions=webdownload&options=offline). It writes the results in VTK format, which can be directly and easily visualized by [Paraview](https://www.paraview.org/).

# Statement of need

Many problems of fundamental importance in solid mechanics involve pre-existing and propagating discontinuities such as cracks. As classical continuum theory employs spatial derivatives in its formulation and assumes the material is continuous as it deforms, it is inherently difficult to predict discontinuous behaviors. In light of the inadequacies of the classical continuum theory, Peridynamics (PD), which is based on non-local interactions and employs spatial integrals, was first introduced in [@silling2000reformulation] for failure analysis of materials and structures.

PD theory has been extended from bond-based theory [@silling2000reformulation], which has a restriction of a fixed Poisson's ratio of 1/4, to the state-based one [@silling_peridynamic_2007], which has no restriction on the value of Poisson's ratio. Both bond-based and state-based PD theories were derived by equating the classical strain energy at a material point to that of PD with a complete sphere interaction domain. However, there PD model requires a surface correction [@le_surface_2018] and a volume correction procedures to improve integration accuracy [@seleson_improved_2014]. To remove these two drawbacks, @madenci_peridynamic_2019 proposed a PD least squares minimization (PDLSM), and @liu_revised_2021 proposed a revised non-ordinary state-based PD. Both methods were derived by Taylor series expansion and the concept of non-local interactions of PD. Comparing to the finite element method (FEM), PD is computationally expensive. Thus, the authors propose a coupling model of PDLSM and FEM (PDLSM-FEM) for taking advantage of these two methods [@liu_simulating_2021;@liu_coupled_2021], and develop the PDLSM-FEM solver in C++ environment.

PDLSM-FEM solver is user friendly and can solve the structures with uniform or non-uniform mesh much faster and save much memory usage comparing to the pure PD model. It does not require surface correction [@le_surface_2018] and volume correction [@seleson_improved_2014] as many published PD models are required. In many published PD theories, constraints are applied through a nonzero volume rather than on a surface, commonly introduced within a fictitious layer [@hu_formulation_2012], and traction at boundaries is introduced in these PD models as body forces within a layer [@madenci_peridynamic_2014]. In the PDLSM-FEM solver, the constraints and traction boundary conditions (BCs) are treated in an easy way as same as FEM, rather than by introducing a fictitious layer. The PDLSM-FEM solver enables 2 criteria for 2D crack propagation simulations: maximum circumferential tensile stress and maximum principal stress.

# Usage and features

PDLSM-FEM is developed to perform fracture analysis of structures, which takes advantage of both PD and FEM. It also can be used to perform conventional finite element analysis, when the PD region of the model is shrunk to zero, which means there is no PD element.
\autoref{fig:diag} illustrates the crack propagation simulation of a diagonal plate with an inclined pre-existing crack by the PDLSM-FEM solver. As it shows, PDLSM-FEM solver captures the crack growth path which has a good agreement with the experimental observation [@ayatollahi_analysis_2009]. The geometry and material properties of this example can be found in [@liu_simulating_2021].

![A diagonal plate under displacement loading. (left-top): Displacement $U_y$ with deformed shape; (right-top): Stress $\sigma_y$ with deformed shape; (left-bottom): Damage $\varphi$ with underformed shape; (right-bottom): Experimental crack path. \label{fig:diag}](diagExample.png){ width=90% }


# References

