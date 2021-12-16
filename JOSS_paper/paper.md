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

PDLSM-FEM solver is a parallel implementation of coupled peridynamics least squares minimization and finite element method (PDLSM-FEM) in 2D and 3D using MPI parallelism. This cross-platform solver is written in a C++ language and includes implicit static, explicit dynamic, and implicit dynamic solvers for structure analysis under displacement or traction loading. PDLSM-FEM solver stores the large sparse matrix in a compressed sparse row (CSR) format rather than a full matrix format, and solves the large linear systems of equations by the Parallel Direct Sparse Solver (PDSS) of [Intel Math Kernel Library](https://software.intel.com/content/www/us/en/develop/tools/oneapi/base-toolkit/download.html?operatingsystem=window&distributions=webdownload&options=offline). It writes the results in VTK format, which can be directly and easily visualized by [Paraview](https://www.paraview.org/).

# Statement of need

Many problems of fundamental importance in solid mechanics involve pre-existing and propagating discontinuities such as cracks. As classical continuum theory employs spatial derivatives in its formulation and assumes the material is continuous as it deforms, it is inherently difficult to predict discontinuous behaviors. In light of the inadequacies of the classical continuum theory, Peridynamics (PD), which is based on non-local interactions and employs spatial integrals, was first introduced in [@silling2000reformulation] for failure analysis of materials and structures.

PD theory has been extended from bond-based theory [@silling2000reformulation], which has a restriction of a fixed Poisson's ratio of 1/4, to the state-based one [@silling_peridynamic_2007], which has no restriction on the value of Poisson's ratio. Both bond-based and state-based PD theories were derived by equating the classical strain energy at a material point to that of PD with a spherical neighborhood. However, these PD models require a surface correction [@le_surface_2018] and a volume correction procedures to improve integration accuracy [@seleson_improved_2014]. To remove these two drawbacks, @madenci_peridynamic_2019 proposed a PD least squares minimization (PDLSM), and @liu_revised_2021 proposed a revised non-ordinary state-based PD. Both methods were derived by Taylor series expansion and the concept of non-local interactions of PD. Comparing to the finite element method (FEM), PD is computationally expensive. Thus, the authors propose a coupling model of PDLSM and FEM (PDLSM-FEM) for taking advantage of these two methods [@liu_simulating_2021;@liu_coupled_2021], and develop the PDLSM-FEM solver in C++ environment.

PDLSM-FEM solver is user friendly and can solve the structures with uniform or non-uniform mesh relatively fast and using less memory compared to pure PD model implementations. It does not require surface correction [@le_surface_2018] and volume correction [@seleson_improved_2014]. For the bond-based PD [@silling2000reformulation] and the state-based PD [@silling_peridynamic_2007], constraints are applied through a nonzero volume rather than on a surface, commonly introduced within a fictitious layer [@hu_formulation_2012], and traction at boundaries is introduced in these PD models as body forces within a layer [@madenci_peridynamic_2014]. In the PDLSM-FEM solver, the constraints and traction boundary conditions (BCs) are treated in an easy way similar to FEM without introducing fictitious layer. The PDLSM-FEM solver enables two criteria only for 2D crack propagation simulations: maximum circumferential tensile stress and maximum principal stress.

# Usage and features

PDLSM-FEM is developed to perform fracture analysis of structures, which takes advantage of both PD and FEM. It can also be used to perform conventional finite element analysis, when the PD region of the model is shrunk to zero, which means there is no PD element.
\autoref{fig:diag} illustrates the crack propagation simulation of a diagonal plate with an inclined pre-existing crack by the PDLSM-FEM solver. As it shows, PDLSM-FEM solver captures the crack growth path which has a good agreement with the experimental observation [@ayatollahi_analysis_2009]. The geometry and material properties of this example can be found in [@liu_simulating_2021].

![Simulation and experimental results of a diagonal plate: displacement $u_y$ with deformed shape (left-top), stress $\sigma_y$ with deformed shape (right-top), damage $\varphi$ with underformed shape (left-bottom), and experimental crack path (right-bottom). \label{fig:diag}](diagExample.png){ width=90% }


# References

