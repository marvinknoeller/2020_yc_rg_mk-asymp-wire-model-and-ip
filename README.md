## License

Copyright (c) 2023, Yves Capdeboscq, Roland Griesmaier, Marvin Knöller

This software is released under GNU General Public License, Version 3.
The full license text is provided in the file LICENSE included in this repository 
or can be obtained from http://www.gnu.org/licenses/


## Content of this project
This is a guide to generate the figures that have been used in the work

„An asymptotic representation formula for scattering by thin tubular structures and an application in inverse scattering“

By Yves Capdeboscq, Roland Griesmaier and Marvin Knöller.

You find all needed Matlab and Python files in the folders. 

## An overview
- [ ] The folder **farfields** includes the far fields corresponding to the three different center curves for different radii that we consider in our work.
- [ ] The folder **grids** includes the corresponding meshes in a .msh file that may be opened using e.g. GMSH (see https://gmsh.info).
- [ ] The folder **matlab_scripts(2020a)** includes all the Matlab files that are needed to compute the plots. 
Details below.
- [ ] The folder **python_bempp_scripts** includes the python scripts for generating the reference solutions.
This requires the boundary element library bempp (see https://bempp.com).
## Matlab script details
In the folder matlab_scripts one finds:
- [ ] the scripts _allpoints.m_, _splineeval.m_, _splinepoints.m_. These are routines for computing with 3d cubic splines.
- [ ] the scripts _Curv_Pen.m_, _Seg_Pen_New.m_, _Jacobian_Curv_Pen.m_, _Jacobian_Seg_Pen_New.m_. These scripts are implementations for the penalty terms.
- [ ] the scripts _FarField_Per_Maxwell_E_spline.m_, _derive_farfield.m_, _SetupFarField.m_, _err_on_ff.m_. Here, we implement the leading order term of the asymptotic perturbation formula and derive it with respect to the center curve. _err_on_ff.m_ computes the (abs. ) error of the far field.
- [ ] the auxiliary scripts _dotReal.m_ and _Tensor_Mat_Mult.m_ implementing the dot product and the matrix multiplication along the third dimension.
- [ ] the plot scripts _Plot_Segment_Comparison.m_ and _Plot_For_Pictures_, plotting Fig.3 and Fig.4 from the paper. These scripts call the _*_Conv_Subseg_paper.m_ and _for*_ scripts.
- [ ] _Reconstruct_examples.m_. This script runs the reconstruction from scratch given the far field from the bempp simulation and eventually plots the results. In this script, change the parameter for the example that you want to reconstruct. This calls the script _inv_test_maxwell_spline.m_.

## Requirements and additional information
The computations have been carried out on a MacBook Pro 15“ model from 2018.
Generating all figures from scratch takes a few hours.
Computations have been carried out using the Matlab 2020a version.

Simulation requires the Matlab parallelization toolbox. 
