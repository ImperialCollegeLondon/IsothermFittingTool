%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Imperial College London, United Kingdom
% Multifunctional Nanomaterials Laboratory / Complex Porous Media
% Laboratory
%
% Project:  PhD
% Year:     2021
% MATLAB:   R2020a
% Authors:  Hassan Azzan (HA)
%
% Purpose:
% Fits input experimental adsorption isotherm data to different isotherm,
% and outputs isotherm parameters and independent confidence bounds
%
% Last modified:
% - 2022-08-08, HA: Add capability to use temperature dependent
%                   statistical model for Zeolites (Ruthven 1982)
% - 2021-11-15, HA: Update confidence region calculation for MLE
% - 2021-11-15, HA: Add capability to use temperature dependent
%                   triple site Langmuir (9 param) isotherm model
% - 2021-11-12, HA: Add capability to use temperature dependent
%                   Dual (9 param) or Single (6 param) site Henry-Langmuir
%                   isotherm model
% - 2021-11-11, HA: Add capability to use 4-parameter temperature dependent
%                   Toth isotherm model
% - 2021-11-10, HA: Add capability to use 6-parameter temperature dependent
%                   virial isotherm model
% - 2021-07-14, HA: Normalise experimetal adsorption data for objective
%                   function calculation
% - 2021-07-08, HA: Add capability save resulting data to a *.mat file
% - 2021-07-01, HA: Add capability to set fixed saturation capacities for
%                   fitting
% - 2021-06-16, HA: Add capability to display isotherm parameters in
%                   concentration units
% - 2021-04-23, HA: Add display of predicted parameters with uncertainty
% - 2021-04-13, HA: Add capability to use single site isotherm models
% - 2021-03-11, HA: Initial creation
%
% Input arguments:
%
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALISATION and INPUTS
% Clear command window and workspace
clc; clear; close all;
% For new data create a new variable (HOME -> New Variable) named 'fitData'
% (rename on workspace) and save the data as *.mat file in 3 column format
% with Pressure (bar), adsorbed amount (-), temperature (K) respectively
% with any name of your choice.
% RUN THIS SCRIPT from ERASE OR IsothermFittingTool folder ONLY!
% Load input experimental data from *.mat or *.csv file via prompt
uiopen
%% Select isotherm model for fitting
% TSL = Triple site Langmuir. DSL = Dual site Langmuir.
% SSL = Single site Langmuir. DSS = Dual site. Sips. SSS = Single site Sips.
% TOTH = Toth Isotherm. TOTH2 = Toth Isotherm (Temp dependent tau). TOTH3 =
% Toth Isotherm (Temp dependent qsat and tau). TOTHCHEM = Toth with
% chemisorption.
% VIRIAL = Virial Equation (4a, 2b). VIRIAL2 = Virial Equation (4a, 4b)
% Henry-DSL = HDSL. Henry-SSL = HSSL, Statistical Zeolite Model = STATZ,
% Statistical Zeolite Model for flexible frameworks= STATZGATE
isothermModel = 'TOTHCHEM';
%% Flag for fitting parameters in concentration units (NOT for virial)
flagConcUnits = 0;
%% Flag for saving output in a matfile (IF TRUE, ENTER FILENAME WHEN PROMPTED IN COMMAND WINDOW)
saveFlag = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       INPUTS COMPLETE. IGNORE THE REST OF THE CODE.                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%            dP             `                                              oo   dP
%            88                                                                 88       
%   88d888b. 88 .d8888b. .d8888b. .d8888b. .d8888b.    dP  dP  dP .d8888b. dP d8888P
%   88'  `88 88 88ooood8 88'  `88 Y8ooooo. 88ooood8    88  88  88 88'  `88 88   88
%   88.  .88 88 88.  ... 88.  .88       88 88.  ...    88.88b.88' 88.  .88 88   88
%   88Y888P' dP `88888P' `88888P8 `88888P' `88888P'    8888P Y8P  `88888P8 dP   dP
%   88
%   dP
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         isotherm parameters can be found in 'isothermData.isothermParameters'          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
isothermData = isothermFittingFun(isothermModel, flagConcUnits, saveFlag, fitData);
toc