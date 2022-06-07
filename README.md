# Adsorption Isotherm fitting tool

Fits input experimental adsorption isotherm data to different isotherm models, and outputs isotherm parameters and parameter confidence intervals.

## How to use
- Open 'isothermFitting.m' and follow instructions at the top of the file.

## Main features

- Fit adsorption isotherms to dual/single-site Langmuir (DSL/SSL), dual/single-site Sips (DSS/SSS), Toth (TOTH), Henry-(dual/single) Langmuir (HDSL/HSSL) and Virial (VIRIAL) isotherms
- Isotherm fitting can be done using either maximum log-likelihood estimation (MLE) or weight least sum of squares (WSS)
- Determines individual confidence intervals based on the diagonal covariance matrix of the of the model to fitted parameters
- Generates contour plots for combinations of parameter pairs to identify correlation between parameters (for Sips and Langmuir models)
- Generates plot of uncertainty in isotherm (solid loading) with respect to pressure based on uncertainty in fitted parameters
- Calculate heat of adsorption as a function of loading using the Virial isotherm fit
- Save fitted parameters and confidence regions in a .mat file

### Other details
- Institution: Imperial College London, United Kingdom
- Multifunctional Nanomaterials Laboratory / Complex Porous Media Laboratory
- Main file: isothermFitting.m
- Project: PhD
- Year: 2021
- MATLAB version: R2020a
- Authors: Hassan Azzan (HA)
- email: hassan.azzan15@imperial.ac.uk

<<<<<<< HEAD
Last updated: 12/08/2021
=======
Last updated: 19/01/2022

## LICENSE
IsothermFittingTool
Copyright (C) 2022  Hassan Azzan

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
>>>>>>> develop
