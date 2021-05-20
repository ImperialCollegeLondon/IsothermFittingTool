# Adsorption Isotherm fitting tool

Fits input experimental adsorption isotherm data to different isotherm models, and outputs isotherm parameters and parameter confidence intervals.

## How to use
- Open 'isothermFitting.m' and follow remaining instructions at the top of the file.

## Main features

- Fit adsorption isotherms to dual/single-site Langmuir (DSL/SSL) and dual/single-site Sips (DSS/SSS) isotherms
- Isotherm fitting can be done using either maximum log-likelihood estimation (MLE) or weight least sum of squares (WSS)
- Determines individual confidence intervals based on the diagonal covariance matrix of the of the model to fitted parameters
- Generates contour plots for combinations of parameter pairs to identify correlation between parameters
- Generates plot of uncertainty in isotherm (solid loading) with respect to pressure based on uncertainty in fitted parameters

### Other details
- Institution: Imperial College London, United Kingdom
- Multifunctional Nanomaterials Laboratory / Complex Porous Media Laboratory
- Main file: isothermFitting.m
- Project: PhD
- Year: 2021
- MATLAB version: R2020a
- Authors: Hassan Azzan (HA)
- email: hassan.azzan15@imperial.ac.uk

Last updated: 20/05/2021
