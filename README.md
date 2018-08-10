# Mullidisciplinary-Structural-Optimization
FEA Analysis of 2D warren Wing Structure with Beam elements 
This Analysis is based on custom FEA analysis of 2D warren truss wing structure using Beam element.
The Structure is modelled start from node points, elements, with Dof using FEA Analysis. 
The material used is Al2024 T3, with density 2780Kg/m^3, Exx 73 GPa, poisons raio 0.33,  analytical static FEA analysis and Modal analysis is done to get the global Stiffness matrix, mass matrix, the deflections, stresses, and frequency values. 
Total 17 elements with Area as a design variable, with 7 constraints are chosen for the sensitivity analysis get the sensitivity data for the given constraints with respect to each design variable. 
The single point function approximation methods are used like linear, reciprocal methods to get the next design point. 
The weight function, constraints, gradients on the concentrated elements are compared.
