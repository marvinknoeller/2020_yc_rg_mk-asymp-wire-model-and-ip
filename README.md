{\rtf1\ansi\ansicpg1252\cocoartf2639
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fmodern\fcharset0 Courier;\f1\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;}
\paperw11900\paperh16840\margl1440\margr1440\vieww15900\viewh15440\viewkind0
\deftab720
\pard\pardeftab720\partightenfactor0

\f0\fs26 \cf0 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 ## License\
\
Copyright (c) 2023, Yves Capdeboscq, Roland Griesmaier, Marvin Kn\'f6ller\
\
This software is released under GNU General Public License, Version 3.\
The full license text is provided in the file LICENSE included in this repository \
or can be obtained from http://www.gnu.org/licenses/\
\
\
## Content of this project\
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0
\cf0 \kerning1\expnd0\expndtw0 \outl0\strokewidth0 This is a guide to generate the figures that have been used in the work\
\
\'84An asymptotic representation formula for scattering by thin tubular structures and an application in inverse scattering\'93\
\
By Yves Capdeboscq, Roland Griesmaier and Marvin Kn\'f6ller.\
\
You find all needed Matlab and Python files in the folders. \
\
## An overview\
\pard\pardeftab720\partightenfactor0
\cf0 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 - [ ]\kerning1\expnd0\expndtw0 \outl0\strokewidth0  The folder farfields includes the far fields corresponding to the three different center curves for different radii that we consider in our work.\
\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 - [ ]\kerning1\expnd0\expndtw0 \outl0\strokewidth0  The folder grids includes the corresponding meshes in a .msh file that may be opened using e.g. GMSH.\
\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 - [ ]\kerning1\expnd0\expndtw0 \outl0\strokewidth0  The folder matlab_scripts(2020a) includes all the Matlab files that are needed to compute the plots. \
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0
\cf0 Details below.\
\pard\pardeftab720\partightenfactor0
\cf0 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 - [ ]\kerning1\expnd0\expndtw0 \outl0\strokewidth0  The folder python_bempp_scripts includes the python scripts for generating the reference solutions.\
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0
\cf0 This requires the boundary element library bempp (see: https://bempp.com ).\
## Matlab script details\
In the folder matlab_scripts one finds:\
\pard\pardeftab720\partightenfactor0
\cf0 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 - [ ]\kerning1\expnd0\expndtw0 \outl0\strokewidth0  the scripts allpoints.m, splineeval.m, splinepoints.m . These are routines for computing with 3d cubic splines.\
\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 - [ ]\kerning1\expnd0\expndtw0 \outl0\strokewidth0  the scripts Curv_Pen.m, Seg_Pen_New.m, Jacobian_Curv_Pen.m, Jacobian_Seg_Pen_New.m . These scripts are implementations for the penalty terms.\
\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 - [ ]\kerning1\expnd0\expndtw0 \outl0\strokewidth0  the scripts FarField_Per_Maxwell_E_spline.m, derive_farfield.m, SetupFarField.m, err_on_ff.m . Here, we implement the leading order term of the asymptotic perturbation formula and derive it with respect to the center curve. err_on_ff.m computes the (abs. ) error of the far field.\
\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 - [ ]\kerning1\expnd0\expndtw0 \outl0\strokewidth0  the auxiliary scripts dotReal.m and Tensor_Mat_Mult.m implementing the dot product and the matrix multiplication along the third dimension.\
\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 - [ ]\kerning1\expnd0\expndtw0 \outl0\strokewidth0  the plot scripts Plot_Segment_Comparison.m and Plot_For_Pictures, plotting Fig.3 and Fig.4 from the paper. These scripts call the *_Conv_Subseg_paper.m and for* scripts.\
\expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 - [ ]\kerning1\expnd0\expndtw0 \outl0\strokewidth0  Reconstruct_examples.m . This script runs the reconstruction from scratch given the far field from the bempp simulation and eventually plots the results. In this script, change the parameter for the example that you want to reconstruct. This calls the script inv_test_maxwell_spline.m .\
\
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0
\cf0 ## Requirements and additional information\
The computations have been carried out on a MacBook Pro 15\'93 model from 2018.\
Generating all figures from scratch takes a few hours.\
Computations have been carried out using the Matlab 2020a version.\
\
Simulation requires the Matlab parallelization toolbox. \

\f1 \
\
}