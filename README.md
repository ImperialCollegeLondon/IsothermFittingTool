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

Last updated: 12/11/2021
