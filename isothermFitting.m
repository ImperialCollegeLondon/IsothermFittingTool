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
% and outputs isotherm parameters and ellipsoidal confidence bounds
%
% Last modified:
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
clc; clear all; close all;
% For new data create a new variable (HOME -> New Variable) named 'fitData'
% (rename on workspace) and save the data as *.mat file in 3 column format
% with Pressure (bar), adsorbed amount (-), temperature (K) respectively
% with any name of your choice.
% RUN THIS SCRIPT from ERASE OR IsothermFittingTool ONLY!!!!!!!!!!
currentDir = strsplit(cd,filesep);
if strcmp(currentDir(end),'ERASE')
    cd IsothermFittingTool
end
% Obtain git commit ID if submodule of ERASE
try
    gitCommitID.ERASE = getGitCommit('..');
    gitCommitID.isotherm = getGitCommit;
    % Else find git commit ID corresponding to isothermFittingTool only
catch
    % Get short version of git commit ID
    [status,cmdout] = system('git rev-parse HEAD');
    if status == 0
        % Command was successful
        gitCommitID.isotherm = cmdout(1:7);
    else
        gitCommitID.isotherm = [];
    end
end
% Load input experimental data from *.mat or *.csv file via prompt
uiopen
% Determine number of bins you want the experimental data to be binned to
% in terms of the total pressure range
% (for Weighted sum of squares method ONLY).
% IF ERROR --> Reduce number of bins until error is gone
nbins = 1;
% Select isotherm model for fitting
% DSL = Dual site Langmuir. SSL = Single site Langmuir. DSS = Dual site
% Sips. SSS = Single site Sips. VIRIAL = Virial Equation
isothermModel = 'VIRIAL';
% Select fitting method.
% WSS = weighted sum of squares, MLE = maximum log-likelihood estimator
% MLE is preferred for data that is from a single source where the error is
% likely to be normally distributed with a mean of 0.
% WSS is preferred for fitting data for cases where the error might not be
% random and not be normally distributed (data from different sources)
fittingMethod = 'MLE';
% Flag for concentration units in parameters
flagConcUnits = 0;
% Flag for plotting objective function contour plots for dual site models
flagContour = 0;
% Flag for plotting statistical plots (q-q plot and error distribution)
flagStats = 0;
% Flag for fixing saturation capacities (0 for CO2 fitting, 1 for other
% gases)
flagFixQsat = 0;
% Flag for saving output in a matfile (IF TRUE, ENTER FILENAME WHEN
% PROMPTED IN COMMAND WINDOW)
saveFlag = 0;
% IF YOU ARE FIXING SATURATION CAPACITY ENTER THE CO2 SATURATION CAPACITIES
% FOR THE RELEVANT MODEL BELOW
qs1 = 1.9834e+00;
qs2 = 8.9856e+00;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       INPUTS COMPLETE. IGNORE THE REST OF THE CODE.                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%            dP                                                            oo   dP
%            88                                                                 88
%   88d888b. 88 .d8888b. .d8888b. .d8888b. .d8888b.    dP  dP  dP .d8888b. dP d8888P
%   88'  `88 88 88ooood8 88'  `88 Y8ooooo. 88ooood8    88  88  88 88'  `88 88   88
%   88.  .88 88 88.  ... 88.  .88       88 88.  ...    88.88b.88' 88.  .88 88   88
%   88Y888P' dP `88888P' `88888P8 `88888P' `88888P'    8888P Y8P  `88888P8 dP   dP
%   88
%   dP
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  isotherm parameters can be found in 'parsDisp' and uncertainties in 'conRange95Disp'  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GENERATING AND SOLVING GLOBAL OPTIMISATION PROBLEM
% Pressure (x), adsorbed amount (z), Temperature (y) data from input
x = fitData(:,1);
z = fitData(:,2);
y = fitData(:,3);
% Reference isotherm parameters for non-dimensionalisation [qs1 qs2 b01 b02 delU1 delU2]
refValsP = [15,15,1e-1,1e-1,4e4,4e4];
refValsC = [15,15,1e-5,1e-5,4e4,4e4];
switch isothermModel
    case 'DSL'
        % Reference isotherm parameters for non-dimensionalisation
        if flagConcUnits
            isoRef = refValsC;
        else
            isoRef = refValsP;
        end
    case 'SSL'
        % Reference isotherm parameters for non-dimensionalisation
        if flagConcUnits
            isoRef = refValsC;
        else
            isoRef = refValsP;
        end
    case  'DSS'
        % Reference isotherm parameters for non-dimensionalisation
        if flagConcUnits
            isoRef = [refValsC 2];
        else
            isoRef = [refValsP 2];
        end
    case  'SSS'
        % Reference isotherm parameters for non-dimensionalisation
        if flagConcUnits
            isoRef = [refValsC 2];
        else
            isoRef = [refValsP 2];
        end
end
% for concentration units, convert pressure to concentration
if ~flagFixQsat
    if flagConcUnits
        x = 1e5*x./(8.314.*y);
    end
    rng default % For reproducibility
    % Create gs, a GlobalSearch solver with its properties set to the defaults.
    gs = GlobalSearch('NumTrialPoints',2400,'NumStageOnePoints',2000,'Display','off');
    % Set fitting procedure based on isotherm model
    switch isothermModel
        case 'DSL'
            % Set objective function based on fitting method
            switch fittingMethod
                case 'MLE'
                    % Number of bins is automatically set to 1 for MLE as
                    % weights cannot be assigned in MLE
                    nbins =1;
                    if length(unique(y)) == 1 % if only one temperature
                        % Generate objective function for MLE method
                        optfunc = @(par) generateMLEfun(x, y, z, nbins, 'DSL', isoRef, par(1), par(2), par(3), ...
                            par(4), 0, 0);
                    else
                        % Generate objective function for MLE method
                        optfunc = @(par) generateMLEfun(x, y, z, nbins, 'DSL', isoRef, par(1), par(2), par(3), ...
                            par(4), par(5), par(6));
                    end
                case 'WSS'
                    % Generate objective function for WSS method
                    if length(unique(y)) == 1
                        % Generate objective function for MLE method
                        optfunc = @(par) generateWSSfun(x, y, z, nbins, 'DSL', isoRef, par(1), par(2), par(3), ...
                            par(4), 0, 0);
                    else
                        % Generate objective function for MLE method
                        optfunc = @(par) generateWSSfun(x, y, z, nbins, 'DSL', isoRef, par(1), par(2), par(3), ...
                            par(4), par(5), par(6));
                    end
            end
            % Initial conditions, lower bounds, and upper bounds for parameters
            % in DSL isotherm model
            if length(unique(y)) == 1 % if only one temperature
                x0 = [0.5, 0.5,0.5,0.5];
                lb = [0,0,0,0];
                ub = [1,1,1,1];
            else
                x0 = [0.5,0.5,0.5,0.5,0.5,0.5];
                lb = [0,0,0,0,0,0];
                ub = [1,1,1,1,1,1];
            end
            % Create global optimisation problem with solver 'fmincon' and
            % other bounds
            problem = createOptimProblem('fmincon','x0',x0,'objective',optfunc,'lb',lb,'ub',ub);
            % Solve the optimisation problem to obtain the isotherm parameters
            % for the fit
            [parVals, fval]= run(gs,problem);
            % Set fitted parameter values for isotherm model calculation
            qs1   = parVals(1).*isoRef(1);
            qs2   = parVals(2).*isoRef(2);
            b01   = parVals(3).*isoRef(3);
            b02   = parVals(4).*isoRef(4);
            if length(unique(y)) == 1 % if only one temperature
                delU1 = 0;
                delU2 = 0;
            else
                delU1 = parVals(5).*isoRef(5);
                delU2 = parVals(6).*isoRef(6);
            end
            % Calculate fitted isotherm loadings for conditions (P,T)
            % corresponding to experimental data
            qfit  = qs1.*(b01.*x.*exp(delU1./(8.314.*y)))./(1+(b01.*x.*exp(delU1./(8.314.*y)))) ...
                + qs2.*(b02.*x.*exp(delU2./(8.314.*y)))./(1+(b02.*x.*exp(delU2./(8.314.*y))));
            % Calculate ellipsoidal confidence intervals (delta parameter) for
            % fitted parameters
            parameters = [qs1, qs2, b01, b02, delU1, delU2];
            parameters(isnan(parameters))=0;
            [conRange95] = conrangeEllipse(x, y, z, qfit, 'DSL', qs1, qs2, b01, b02, delU1, delU2);
            conRange95(isnan(conRange95))=0;
            % Convert confidence intervals to percentage error
            percentageError = conRange95./parameters' *100;
            fprintf('Isotherm model: %s \n', isothermModel);
            if ~flagConcUnits
                parNames = ["qs1" "qs2" "b01" "b02" "delU1" "delU2"];
                units = ["mol/kg" "mol/kg" "1/bar" "1/bar" "J/mol" "J/mol"];
                parsDisp = parameters;
                conRange95Disp = conRange95;
                for ii = 1:length(parameters)
                    if parsDisp(ii) == 0
                    else
                        fprintf('%s = %5.4e ± %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            else
                parNames = ["qs1" "qs2" "b01" "b02" "delU1" "delU2"];
                units = ["mol/kg" "mol/kg" "m3/mol" "m3/mol" "J/mol" "J/mol"];
                parsDisp = parameters;
                conRange95Disp = conRange95;
                for ii = 1:length(parsDisp)
                    if parsDisp(ii) == 0
                    else
                        fprintf('%s = %5.4e ± %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            end
        case 'SSL'
            % Set objective function based on fitting method
            switch fittingMethod
                case 'MLE'
                    % Number of bins is automatically set to 1 for MLE as
                    % weights cannot be assigned in MLE
                    nbins =1;
                    if length(unique(y)) == 1 % if only one temperature
                        % Generate objective function for MLE method
                        optfunc = @(par) generateMLEfun(x, y, z, nbins, 'DSL', isoRef, par(1), 0, par(2), ...
                            0, 0, 0);
                    else
                        % Generate objective function for MLE method
                        optfunc = @(par) generateMLEfun(x, y, z, nbins, 'DSL', isoRef, par(1), 0, par(2), ...
                            0, par(3), 0);
                    end
                case 'WSS'
                    % Generate objective function for WSS method
                    if length(unique(y)) == 1
                        % Generate objective function for MLE method
                        optfunc = @(par) generateWSSfun(x, y, z, nbins, 'DSL', isoRef, par(1), 0, par(2), ...
                            0, 0, 0);
                    else
                        % Generate objective function for MLE method
                        optfunc = @(par) generateWSSfun(x, y, z, nbins, 'DSL', isoRef, par(1), 0, par(2), ...
                            0, par(3), 0);
                    end
            end
            % Initial conditions, lower bounds, and upper bounds for parameters
            % in DSL isotherm model
            if length(unique(y)) == 1 % if only one temperature
                x0 = [0.5,0.5];
                lb = [0,0];
                ub = [1,1];
            else
                x0 = [0.5,0.5,0.5];
                lb = [0,0,0];
                ub = [1,1,1];
            end
            % Create global optimisation problem with solver 'fmincon' and
            % other bounds
            problem = createOptimProblem('fmincon','x0',x0,'objective',optfunc,'lb',lb,'ub',ub);
            % Solve the optimisation problem to obtain the isotherm parameters
            % for the fit
            [parVals, fval]= run(gs,problem);
            % Set fitted parameter values for isotherm model calculation
            qs1   = parVals(1).*isoRef(1);
            qs2   = 0;
            b01   = parVals(2).*isoRef(3);
            b02   = 0;
            if length(unique(y)) == 1 % if only one temperature
                delU1 = 0;
            else
                delU1 = parVals(3).*isoRef(5);
            end
            delU2 = 0;
            % Calculate fitted isotherm loadings for conditions (P,T)
            % corresponding to experimental data
            qfit  = qs1.*(b01.*x.*exp(delU1./(8.314.*y)))./(1+(b01.*x.*exp(delU1./(8.314.*y)))) ...
                + qs2.*(b02.*x.*exp(delU2./(8.314.*y)))./(1+(b02.*x.*exp(delU2./(8.314.*y))));
            % Calculate ellipsoidal confidence intervals (delta parameter) for
            % fitted parameters
            parameters = [qs1, qs2, b01, b02, delU1, delU2];
            parameters(isnan(parameters))=0;
            [conRange95] = conrangeEllipse(x, y, z, qfit, 'DSL', qs1, qs2, b01, b02, delU1, delU2);
            conRange95(isnan(conRange95))=0;
            % Convert confidence intervals to percentage error
            percentageError = conRange95./parameters' *100;
            fprintf('Isotherm model: %s \n', isothermModel);
            if ~flagConcUnits
                parNames = ["qs1" "qs2" "b01" "b02" "delU1" "delU2"];
                units = ["mol/kg" "mol/kg" "1/bar" "1/bar" "J/mol" "J/mol"];
                parsDisp = parameters;
                conRange95Disp = conRange95;
                for ii = 1:length(parameters)
                    if parsDisp(ii) == 0
                    else
                        fprintf('%s = %5.4e ± %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            else
                parNames = ["qs1" "qs2" "b01" "b02" "delU1" "delU2"];
                units = ["mol/kg" "mol/kg" "m3/mol" "m3/mol" "J/mol" "J/mol"];
                parsDisp = parameters;
                conRange95Disp = conRange95;
                for ii = 1:length(parsDisp)
                    if parsDisp(ii) == 0
                    else
                        fprintf('%s = %5.4e ± %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            end
        case 'DSS'
            % Set objective function based on fitting method
            switch fittingMethod
                case 'MLE'
                    % Number of bins is automatically set to 1 for MLE as
                    % weights cannot be assigned in MLE
                    nbins =1;
                    % Generate objective function for MLE method
                    optfunc = @(par) generateMLEfun(x, y, z, nbins, isothermModel, isoRef, par(1), par(2), par(3), ...
                        par(4), par(5), par(6), par(7));
                case 'WSS'
                    % Generate objective function for WSS method
                    optfunc = @(par) generateWSSfun(x, y, z, nbins, isothermModel, isoRef, par(1), par(2), par(3), ...
                        par(4), par(5), par(6), par(7));
            end
            % Initial conditions, lower bounds, and upper bounds for parameters
            % in DSL isotherm model
            x0 = [0.5,0.5,0.5,0.5,0.5,0.5,0.5];
            lb = [0,0,0,0,0,0,0];
            ub = [1,1,1,1,1,1,1];
            % Create global optimisation problem with solver 'fmincon' and
            % other bounds
            problem = createOptimProblem('fmincon','x0',x0,'objective',optfunc,'lb',lb,'ub',ub);
            % Solve the optimisation problem to obtain the isotherm parameters
            % for the fit
            [parVals, fval]= run(gs,problem);
            % Set fitted parameter values for isotherm model calculation
            qs1   = parVals(1).*isoRef(1);
            qs2   = parVals(2).*isoRef(2);
            b01   = parVals(3).*isoRef(3);
            b02   = parVals(4).*isoRef(4);
            delU1 = parVals(5).*isoRef(5);
            delU2 = parVals(6).*isoRef(6);
            gamma = parVals(7).*isoRef(7);
            % Calculate fitted isotherm loadings for conditions (P,T)
            % corresponding to experimental data
            qfit  = qs1.*(b01.*x.*exp(delU1./(8.314.*y))).^gamma./(1+(b01.*x.*exp(delU1./(8.314.*y))).^gamma) ...
                + qs2.*(b02.*x.*exp(delU2./(8.314.*y))).^gamma./(1+(b02.*x.*exp(delU2./(8.314.*y))).^gamma);
            % Calculate ellipsoidal confidence intervals (delta parameter) for
            % fitted parameters
            parameters = [qs1, qs2, b01, b02, delU1, delU2, gamma];
            parameters(isnan(parameters))=0;
            [conRange95] = conrangeEllipse(x, y, z, qfit, isothermModel, qs1, qs2, b01, b02, delU1, delU2, gamma);
            conRange95(isnan(conRange95))=0;
            % Convert confidence intervals to percentage error
            percentageError = conRange95./parameters' *100;
            fprintf('Isotherm model: %s \n', isothermModel);
            if ~flagConcUnits
                parNames = ["qs1" "qs2" "b01" "b02" "delU1" "delU2" "gamma"];
                units = ["mol/kg" "mol/kg" "1/bar" "1/bar" "J/mol" "J/mol" " "];
                parsDisp = parameters;
                conRange95Disp = conRange95;
                for ii = 1:length(parameters)
                    if parsDisp(ii) == 0
                    else
                        fprintf('%s = %5.4e ± %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            else
                parNames = ["qs1" "qs2" "b01" "b02" "delU1" "delU2"  "gamma"];
                units = ["mol/kg" "mol/kg" "m3/mol" "m3/mol" "J/mol" "J/mol"  " "];
                parsDisp = parameters;
                conRange95Disp = conRange95;
                for ii = 1:length(parsDisp)
                    if parsDisp(ii) == 0
                    else
                        fprintf('%s = %5.4e ± %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            end
        case 'SSS'
            % Set objective function based on fitting method
            switch fittingMethod
                case 'MLE'
                    % Number of bins is automatically set to 1 for MLE as
                    % weights cannot be assigned in MLE
                    nbins =1;
                    % Generate objective function for MLE method
                    optfunc = @(par) generateMLEfun(x, y, z, nbins, 'DSS', isoRef, par(1), 0, par(2), ...
                        0, par(3), 0, par(4));
                case 'WSS'
                    % Generate objective function for WSS method
                    optfunc = @(par) generateWSSfun(x, y, z, nbins, 'DSS', isoRef, par(1), 0, par(2), ...
                        0, par(3), 0, par(4));
            end
            % Initial conditions, lower bounds, and upper bounds for parameters
            % in DSL isotherm model
            x0 = [0.5,0.5,0.5,0.5];
            lb = [0,0,0,0];
            ub = [1,1,1,1];
            % Create global optimisation problem with solver 'fmincon' and
            % other bounds
            problem = createOptimProblem('fmincon','x0',x0,'objective',optfunc,'lb',lb,'ub',ub);
            % Solve the optimisation problem to obtain the isotherm parameters
            % for the fit
            [parVals, fval]= run(gs,problem);
            % Set fitted parameter values for isotherm model calculation
            qs1   = parVals(1).*isoRef(1);
            qs2   = 0;
            b01   = parVals(2).*isoRef(3);
            b02   = 0;
            delU1 = parVals(3).*isoRef(5);
            delU2 = 0;
            gamma = parVals(4).*isoRef(7);
            % Calculate fitted isotherm loadings for conditions (P,T)
            % corresponding to experimental data
            qfit  = qs1.*(b01.*x.*exp(delU1./(8.314.*y))).^gamma./(1+(b01.*x.*exp(delU1./(8.314.*y))).^gamma) ...
                + qs2.*(b02.*x.*exp(delU2./(8.314.*y))).^gamma./(1+(b02.*x.*exp(delU2./(8.314.*y))).^gamma);
            % Calculate ellipsoidal confidence intervals (delta parameter) for
            % fitted parameters
            parameters = [qs1, qs2, b01, b02, delU1, delU2, gamma];
            parameters(isnan(parameters))=0;
            [conRange95] = conrangeEllipse(x, y, z, qfit, 'DSS', qs1, qs2, b01, b02, delU1, delU2, gamma);
            conRange95(isnan(conRange95))=0;
            % Convert confidence intervals to percentage error
            percentageError = conRange95./parameters' *100;
            fprintf('Isotherm model: %s \n', isothermModel);
            if ~flagConcUnits
                parNames = ["qs1" "qs2" "b01" "b02" "delU1" "delU2" "gamma"];
                units = ["mol/kg" "mol/kg" "1/bar" "1/bar" "J/mol" "J/mol" " "];
                parsDisp = parameters;
                conRange95Disp = conRange95;
                for ii = 1:length(parameters)
                    if parsDisp(ii) == 0
                    else
                        fprintf('%s = %5.4e ± %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            else
                parNames = ["qs1" "qs2" "b01" "b02" "delU1" "delU2"  "gamma"];
                units = ["mol/kg" "mol/kg" "m3/mol" "m3/mol" "J/mol" "J/mol"  " "];
                parsDisp = parameters;
                conRange95Disp = conRange95;
                for ii = 1:length(parsDisp)
                    if parsDisp(ii) == 0
                    else
                        fprintf('%s = %5.4e ± %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            end
        case 'VIRIAL'
            % Reference isotherm parameters for non-dimensionalisation [qs1 qs2 b01 b02 delU1 delU2]
            refValsP = [10e4,10e4,10e4,10e4,10e4,10e4,10e4,10e4];
            refValsC = refValsP;
            isoRef = refValsC;
            % Set objective function based on fitting method
            switch fittingMethod
                case 'MLE'
                    % Number of bins is automatically set to 1 for MLE as
                    % weights cannot be assigned in MLE
                    nbins =1;
                    % Generate objective function for MLE method
                    optfunc = @(par) generateMLEfun(x, y, z, nbins, 'VIRIAL', isoRef, par(1), par(2), ...
                        par(3), par(4), par(5), par(6), 0, 0);
                case 'WSS'
                    % Generate objective function for WSS method
                    optfunc = @(par) generateWSSfun(x, y, z, nbins, 'VIRIAL', isoRef, par(1), par(2), ...
                        par(3), par(4), par(5), par(6), 0, 0);
            end
            % Initial conditions, lower bounds, and upper bounds for parameters
            % in DSL isotherm model
            x0 = [0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
            lb = [-1,-1,-1,-1,-1,-1,-1,-1];
            ub = [1,1,1,1,1,1,1,1];
            % Create global optimisation problem with solver 'fmincon' and
            % other bounds
            problem = createOptimProblem('fmincon','x0',x0,'objective',optfunc,'lb',lb,'ub',ub);
            % Solve the optimisation problem to obtain the isotherm parameters
            % for the fit
            [parVals, fval]= run(gs,problem);
            % Set fitted parameter values for isotherm model calculation
            a0   = parVals(1).*isoRef(1);
            a1   = parVals(2).*isoRef(2);
            a2   = parVals(3).*isoRef(3);
            a3   = parVals(4).*isoRef(4);
            b0   = parVals(5).*isoRef(5);
            b1   = parVals(6).*isoRef(6);
            b2   = 0;
            b3   = 0;
            parameters = [a0, a1, a2, a3, b0, b1, b2, b3];
            parameters(isnan(parameters))=0;
            % Calculate fitted isotherm loadings for conditions (P,T)
            % corresponding to experimental data
            qfit = [];
            expData = [x,z,y];
            expData = sortrows(expData,1);
            expData = expData(length(unique(y))+1:end,:);
            
            x = expData(:,1);
            z = expData(:,2);
            y = expData(:,3);
            for kk = 1:length(x)
                lnPfit(kk)  = log(z(kk)) + 1/y(kk).*(parameters(1) + parameters(2)*z(kk) + ...
                    parameters(3)*z(kk)^2 + parameters(4)*z(kk)^3)  + parameters(5) ...
                    + parameters(6)*z(kk)+ parameters(7)*z(kk).^2 ...
                    + parameters(8)*z(kk).^3;
            end
            % Calculate ellipsoidal confidence intervals (delta parameter) for
            % fitted parameters
            [conRange95] = conrangeEllipse(x, y, z, lnPfit, 'VIRIAL', a0, a1, a2, a3, b0, b1, b2, b3);
            conRange95(isnan(conRange95))=0;
            % Convert confidence intervals to percentage error
            percentageError = conRange95./parameters' *100;
            fprintf('Isotherm model: %s \n', isothermModel);
            parNames = ["a0" "a1" "a2" "a3" "b0" "b1" "b2" "b3"];
            units = ["K" "1/K" "1/K^2" "1/K^3" " " " " " " " "];
            parsDisp = parameters;
            conRange95Disp = conRange95;
            for ii = 1:length(parameters)
                if parsDisp(ii) == 0
                else
                    fprintf('%s = %5.4e ± %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                end
            end
    end
else
    if flagConcUnits
        x = 1e5*x./(8.314.*y);
    end
    rng default % For reproducibility
    % Create gs, a GlobalSearch solver with its properties set to the defaults.
    gs = GlobalSearch('NumTrialPoints',1400,'NumStageOnePoints',400,'Display','off');
    % Set fitting procedure based on isotherm model
    switch isothermModel
        case 'DSL'
            % Set objective function based on fitting method
            switch fittingMethod
                case 'MLE'
                    % Number of bins is automatically set to 1 for MLE as
                    % weights cannot be assigned in MLE
                    nbins =1;
                    if length(unique(y)) == 1 % if only one temperature
                        % Generate objective function for MLE method
                        optfunc = @(par) generateMLEfun(x, y, z, nbins, 'DSL', isoRef, qs1, qs2, par(1), ...
                            par(2), 0, 0);
                    else
                        % Generate objective function for MLE method
                        optfunc = @(par) generateMLEfun(x, y, z, nbins, 'DSL', isoRef, qs1, qs2, par(1), ...
                            par(2), par(3), par(4));
                    end
                case 'WSS'
                    % Generate objective function for WSS method
                    if length(unique(y)) == 1
                        % Generate objective function for MLE method
                        optfunc = @(par) generateWSSfun(x, y, z, nbins, 'DSL', isoRef, qs1, qs2, par(1), ...
                            par(2), 0, 0);
                    else
                        % Generate objective function for MLE method
                        optfunc = @(par) generateWSSfun(x, y, z, nbins, 'DSL', isoRef, qs1, qs2, par(1), ...
                            par(2), par(3), par(4));
                    end
            end
            % Initial conditions, lower bounds, and upper bounds for parameters
            % in DSL isotherm model
            if length(unique(y)) == 1 % if only one temperature
                x0 = [0.5,0.5];
                lb = [0,0];
                ub = [1,1];
            else
                x0 = [0.5,0.5,0.5,0.5];
                lb = [0,0,0,0];
                ub = [1,1,1,1];
            end
            % Create global optimisation problem with solver 'fmincon' and
            % other bounds
            problem = createOptimProblem('fmincon','x0',x0,'objective',optfunc,'lb',lb,'ub',ub);
            % Solve the optimisation problem to obtain the isotherm parameters
            % for the fit
            [parVals, fval]= run(gs,problem);
            % Set fitted parameter values for isotherm model calculation
            qs1   = qs1;
            qs2   = qs2;
            b01   = parVals(1).*isoRef(3);
            b02   = parVals(2).*isoRef(4);
            if length(unique(y)) == 1 % if only one temperature
                delU1 = 0;
                delU2 = 0;
            else
                delU1 = parVals(3).*isoRef(5);
                delU2 = parVals(4).*isoRef(6);
            end
            % Calculate fitted isotherm loadings for conditions (P,T)
            % corresponding to experimental data
            qfit  = qs1.*(b01.*x.*exp(delU1./(8.314.*y)))./(1+(b01.*x.*exp(delU1./(8.314.*y)))) ...
                + qs2.*(b02.*x.*exp(delU2./(8.314.*y)))./(1+(b02.*x.*exp(delU2./(8.314.*y))));
            % Calculate ellipsoidal confidence intervals (delta parameter) for
            % fitted parameters
            parameters = [qs1, qs2, b01, b02, delU1, delU2];
            parameters(isnan(parameters))=0;
            [conRange95] = conrangeEllipse(x, y, z, qfit, 'DSL', qs1, qs2, b01, b02, delU1, delU2);
            conRange95(isnan(conRange95))=0;
            % Convert confidence intervals to percentage error
            percentageError = conRange95./parameters' *100;
            if ~flagConcUnits
                parNames = ["qs1" "qs2" "b01" "b02" "delU1" "delU2"];
                units = ["mol/kg" "mol/kg" "1/bar" "1/bar" "J/mol" "J/mol"];
                parsDisp = parameters;
                conRange95Disp = conRange95;
                for ii = 1:length(parameters)
                    if parsDisp(ii) == 0
                    else
                        fprintf('%s = %5.4e ± %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            else
                parNames = ["qs1" "qs2" "b01" "b02" "delU1" "delU2"];
                units = ["mol/kg" "mol/kg" "m3/mol" "m3/mol" "J/mol" "J/mol"];
                parsDisp = parameters;
                conRange95Disp = conRange95;
                for ii = 1:length(parsDisp)
                    if parsDisp(ii) == 0
                    else
                        fprintf('%s = %5.4e ± %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            end
        case 'SSL'
            % Set objective function based on fitting method
            switch fittingMethod
                case 'MLE'
                    % Number of bins is automatically set to 1 for MLE as
                    % weights cannot be assigned in MLE
                    nbins =1;
                    if length(unique(y)) == 1 % if only one temperature
                        % Generate objective function for MLE method
                        optfunc = @(par) generateMLEfun(x, y, z, nbins, 'DSL', isoRef, qs1, 0, par(1), ...
                            0, 0, 0);
                    else
                        % Generate objective function for MLE method
                        optfunc = @(par) generateMLEfun(x, y, z, nbins, 'DSL', isoRef, qs1, 0, par(1), ...
                            0, par(2), 0);
                    end
                case 'WSS'
                    % Generate objective function for WSS method
                    if length(unique(y)) == 1
                        % Generate objective function for MLE method
                        optfunc = @(par) generateWSSfun(x, y, z, nbins, 'DSL', isoRef, qs1, 0, par(1), ...
                            0, 0, 0);
                    else
                        % Generate objective function for MLE method
                        optfunc = @(par) generateWSSfun(x, y, z, nbins, 'DSL', isoRef, qs1, 0, par(1), ...
                            0, par(2), 0);
                    end
            end
            % Initial conditions, lower bounds, and upper bounds for parameters
            % in DSL isotherm model
            if length(unique(y)) == 1 % if only one temperature
                x0 = [0.5];
                lb = [0];
                ub = [1];
            else
                x0 = [0.5,0.5];
                lb = [0,0];
                ub = [1,1];
            end
            % Create global optimisation problem with solver 'fmincon' and
            % other bounds
            problem = createOptimProblem('fmincon','x0',x0,'objective',optfunc,'lb',lb,'ub',ub);
            % Solve the optimisation problem to obtain the isotherm parameters
            % for the fit
            [parVals, fval]= run(gs,problem);
            % Set fitted parameter values for isotherm model calculation
            qs1   = qs1;
            qs2   = 0;
            b01   = parVals(1).*isoRef(3);
            b02   = 0;
            if length(unique(y)) == 1 % if only one temperature
                delU1 = 0;
            else
                delU1 = parVals(2).*isoRef(5);
            end
            delU2 = 0;
            % Calculate fitted isotherm loadings for conditions (P,T)
            % corresponding to experimental data
            qfit  = qs1.*(b01.*x.*exp(delU1./(8.314.*y)))./(1+(b01.*x.*exp(delU1./(8.314.*y)))) ...
                + qs2.*(b02.*x.*exp(delU2./(8.314.*y)))./(1+(b02.*x.*exp(delU2./(8.314.*y))));
            % Calculate ellipsoidal confidence intervals (delta parameter) for
            % fitted parameters
            parameters = [qs1, qs2, b01, b02, delU1, delU2];
            parameters(isnan(parameters))=0;
            [conRange95] = conrangeEllipse(x, y, z, qfit, 'DSL', qs1, qs2, b01, b02, delU1, delU2);
            conRange95(isnan(conRange95))=0;
            % Convert confidence intervals to percentage error
            percentageError = conRange95./parameters' *100;
            if ~flagConcUnits
                parNames = ["qs1" "qs2" "b01" "b02" "delU1" "delU2"];
                units = ["mol/kg" "mol/kg" "1/bar" "1/bar" "J/mol" "J/mol"];
                parsDisp = parameters;
                conRange95Disp = conRange95;
                for ii = 1:length(parameters)
                    if parsDisp(ii) == 0
                    else
                        fprintf('%s = %5.4e ± %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            else
                parNames = ["qs1" "qs2" "b01" "b02" "delU1" "delU2"];
                units = ["mol/kg" "mol/kg" "m3/mol" "m3/mol" "J/mol" "J/mol"];
                parsDisp = parameters;
                conRange95Disp = conRange95;
                for ii = 1:length(parsDisp)
                    if parsDisp(ii) == 0
                    else
                        fprintf('%s = %5.4e ± %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            end
        case 'DSS'
            % Set objective function based on fitting method
            switch fittingMethod
                case 'MLE'
                    % Number of bins is automatically set to 1 for MLE as
                    % weights cannot be assigned in MLE
                    nbins =1;
                    % Generate objective function for MLE method
                    optfunc = @(par) generateMLEfun(x, y, z, nbins, isothermModel, isoRef, qs1, qs2, par(1), ...
                        par(2), par(3), par(4), par(5));
                case 'WSS'
                    % Generate objective function for WSS method
                    optfunc = @(par) generateWSSfun(x, y, z, nbins, isothermModel, isoRef, qs1, qs2, par(1), ...
                        par(2), par(3), par(4), par(5));
            end
            % Initial conditions, lower bounds, and upper bounds for parameters
            % in DSL isotherm model
            x0 = [0.5,0.5,0.5,0.5,0.5];
            lb = [0,0,0,0,0];
            ub = [1,1,1,1,1];
            % Create global optimisation problem with solver 'fmincon' and
            % other bounds
            problem = createOptimProblem('fmincon','x0',x0,'objective',optfunc,'lb',lb,'ub',ub);
            % Solve the optimisation problem to obtain the isotherm parameters
            % for the fit
            [parVals, fval]= run(gs,problem);
            % Set fitted parameter values for isotherm model calculation
            qs1   = qs1;
            qs2   = qs2;
            b01   = parVals(1).*isoRef(3);
            b02   = parVals(2).*isoRef(4);
            delU1 = parVals(3).*isoRef(5);
            delU2 = parVals(4).*isoRef(6);
            gamma = parVals(5).*isoRef(7);
            % Calculate fitted isotherm loadings for conditions (P,T)
            % corresponding to experimental data
            qfit  = qs1.*(b01.*x.*exp(delU1./(8.314.*y))).^gamma./(1+(b01.*x.*exp(delU1./(8.314.*y))).^gamma) ...
                + qs2.*(b02.*x.*exp(delU2./(8.314.*y))).^gamma./(1+(b02.*x.*exp(delU2./(8.314.*y))).^gamma);
            % Calculate ellipsoidal confidence intervals (delta parameter) for
            % fitted parameters
            parameters = [qs1, qs2, b01, b02, delU1, delU2, gamma];
            parameters(isnan(parameters))=0;
            [conRange95] = conrangeEllipse(x, y, z, qfit, isothermModel, qs1, qs2, b01, b02, delU1, delU2, gamma);
            conRange95(isnan(conRange95))=0;
            % Convert confidence intervals to percentage error
            percentageError = conRange95./parameters' *100;
            if ~flagConcUnits
                parNames = ["qs1" "qs2" "b01" "b02" "delU1" "delU2" "gamma"];
                units = ["mol/kg" "mol/kg" "1/bar" "1/bar" "J/mol" "J/mol" " "];
                parsDisp = parameters;
                conRange95Disp = conRange95;
                for ii = 1:length(parameters)
                    if parsDisp(ii) == 0
                    else
                        fprintf('%s = %5.4e ± %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            else
                parNames = ["qs1" "qs2" "b01" "b02" "delU1" "delU2"  "gamma"];
                units = ["mol/kg" "mol/kg" "m3/mol" "m3/mol" "J/mol" "J/mol"  " "];
                parsDisp = parameters;
                conRange95Disp = conRange95;
                for ii = 1:length(parsDisp)
                    if parsDisp(ii) == 0
                    else
                        fprintf('%s = %5.4e ± %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            end
        case 'SSS'
            % Set objective function based on fitting method
            switch fittingMethod
                case 'MLE'
                    % Number of bins is automatically set to 1 for MLE as
                    % weights cannot be assigned in MLE
                    nbins =1;
                    % Generate objective function for MLE method
                    optfunc = @(par) generateMLEfun(x, y, z, nbins, 'DSS', isoRef, qs1, 0, par(1), ...
                        0, par(2), 0, par(3));
                case 'WSS'
                    % Generate objective function for WSS method
                    optfunc = @(par) generateWSSfun(x, y, z, nbins, 'DSS', isoRef, qs1, 0, par(1), ...
                        0, par(2), 0, par(3));
            end
            % Initial conditions, lower bounds, and upper bounds for parameters
            % in DSL isotherm model
            x0 = [0.5,0.5,0.5];
            lb = [0,0,0];
            ub = [1,1,1];
            % Create global optimisation problem with solver 'fmincon' and
            % other bounds
            problem = createOptimProblem('fmincon','x0',x0,'objective',optfunc,'lb',lb,'ub',ub);
            % Solve the optimisation problem to obtain the isotherm parameters
            % for the fit
            [parVals, fval]= run(gs,problem);
            % Set fitted parameter values for isotherm model calculation
            qs1   = qs1;
            qs2   = 0;
            b01   = parVals(1).*isoRef(3);
            b02   = 0;
            delU1 = parVals(2).*isoRef(5);
            delU2 = 0;
            gamma = parVals(3).*isoRef(7);
            % Calculate fitted isotherm loadings for conditions (P,T)
            % corresponding to experimental data
            qfit  = qs1.*(b01.*x.*exp(delU1./(8.314.*y))).^gamma./(1+(b01.*x.*exp(delU1./(8.314.*y))).^gamma) ...
                + qs2.*(b02.*x.*exp(delU2./(8.314.*y))).^gamma./(1+(b02.*x.*exp(delU2./(8.314.*y))).^gamma);
            % Calculate ellipsoidal confidence intervals (delta parameter) for
            % fitted parameters
            parameters = [qs1, qs2, b01, b02, delU1, delU2, gamma];
            parameters(isnan(parameters))=0;
            [conRange95] = conrangeEllipse(x, y, z, qfit, 'DSS', qs1, qs2, b01, b02, delU1, delU2, gamma);
            conRange95(isnan(conRange95))=0;
            % Convert confidence intervals to percentage error
            percentageError = conRange95./parameters' *100;
            if ~flagConcUnits
                parNames = ["qs1" "qs2" "b01" "b02" "delU1" "delU2" "gamma"];
                units = ["mol/kg" "mol/kg" "1/bar" "1/bar" "J/mol" "J/mol" " "];
                parsDisp = parameters;
                conRange95Disp = conRange95;
                for ii = 1:length(parameters)
                    if parsDisp(ii) == 0
                    else
                        fprintf('%s = %5.4e ± %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            else
                parNames = ["qs1" "qs2" "b01" "b02" "delU1" "delU2"  "gamma"];
                units = ["mol/kg" "mol/kg" "m3/mol" "m3/mol" "J/mol" "J/mol"  " "];
                parsDisp = parameters;
                conRange95Disp = conRange95;
                for ii = 1:length(parsDisp)
                    if parsDisp(ii) == 0
                    else
                        fprintf('%s = %5.4e ± %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            end
    end
end
fprintf('%s %5.4e \n','objective function:',fval);
%% PLOT RESULTING OUTPUTS
switch isothermModel
    case 'DSL'
        [outScatter,uncBounds]=generateUncertaintySpread(x,y,z,'DSL',parameters,conRange95);
    case 'SSL'
        [outScatter,uncBounds]=generateUncertaintySpread(x,y,z,'DSL',parameters,conRange95);
    case 'DSS'
        [outScatter,uncBounds]=generateUncertaintySpread(x,y,z,'DSS',parameters,conRange95);
    case 'SSS'
        [outScatter,uncBounds]=generateUncertaintySpread(x,y,z,'DSS',parameters,conRange95);
    case 'VIRIAL'
        [outScatter,uncBounds]=generateUncertaintySpread(x,y,z,'VIRIAL',parameters,conRange95);
end
switch isothermModel
    case 'VIRIAL'
        qvals = linspace(0,max(z),500);
        Tvals = unique(y);
        lnPvals = zeros(length(qvals),length(Tvals));
        for jj = 1:length(qvals)
            for kk = 1:length(Tvals)
                qeq = qvals(jj);
                T = Tvals(kk);
                lnPvals(jj,kk) = log(qeq) + 1/T.*(parameters(1) + parameters(2).*qeq + ...
                    parameters(3).*qeq^2 + parameters(4).*qeq^3) + parameters(5) ...
                    + parameters(6).*qeq+ parameters(7).*qeq.^2 ...
                    + parameters(8).*qeq.^3;
            end
        end
        
        figure(1)
        subplot(1,3,1)
        scatter(uncBounds(:,1),uncBounds(:,2),0.5,'MarkerEdgeColor','b','MarkerEdgeAlpha',0.2)
        hold on
        scatter(uncBounds(:,1),uncBounds(:,3),0.5,'MarkerEdgeColor','b','MarkerEdgeAlpha',0.2)
        for kk = 1:length(Tvals)
            plot(qvals,lnPvals(:,kk),'-k','LineWidth',1.5);
        end
        plot(z,log(x),'ob');
        xlabel('Amount adsorbed [mol/kg]');
        ylabel('ln(P) [-]');
        % quantile-quantile plot of experimental data vs normal distribution
        box on
        set(gca,'YScale','linear','XScale','linear','FontSize',15,'LineWidth',1)
        grid on; axis square
        subplot(1,3,2)
        hold on
        for kk = 1:length(Tvals)
            plot(exp(lnPvals(:,kk)),qvals,'-k','LineWidth',1.5);
        end
        plot(x,z,'ob');
        xlabel('Pressure [bar]');
        ylabel('Amount adsorbed [mol/kg]');
        xlim([0 max(x)+1]);
        ylim([0 max(z)+1]);
        box on
        set(gca,'YScale','linear','XScale','linear','FontSize',15,'LineWidth',1)
        grid on; axis square
        subplot(1,3,3)
        for kk = 1:length(Tvals)
            delH = -8.314.*(a0 +a1.*qvals+a2.*qvals.^2+a3.*qvals.^3)./1000;
            plot(qvals,delH,'-k','LineWidth',1.5);
        end
        xlabel('Amount adsorbed [mol/kg]');
        ylabel('-\Delta H_{ads} (Heat of Adsorption) [kJ/mol]');
        ylim([0 max(delH)+5]);
        box on
        set(gca,'YScale','linear','XScale','linear','FontSize',15,'LineWidth',1)
        grid on; axis square
        
        outFit = [qvals' lnPvals];
    otherwise
        % plot of experimental data and fitted data (q vs P)
        Pvals = linspace(0,max(x),500);
        Tvals = unique(y);
        qvals = zeros(length(Pvals),length(Tvals));
        for jj = 1:length(Pvals)
            for kk = 1:length(Tvals)
                P = Pvals(jj);
                T = Tvals(kk);
                switch isothermModel
                    case 'DSL'
                        qvals(jj,kk) = qs1.*(b01.*P.*exp(delU1./(8.314.*T)))./(1+(b01.*P.*exp(delU1./(8.314.*T)))) ...
                            + qs2.*(b02.*P.*exp(delU2./(8.314.*T)))./(1+(b02.*P.*exp(delU2./(8.314.*T))));
                    case 'SSL'
                        qvals(jj,kk) = qs1.*(b01.*P.*exp(delU1./(8.314.*T)))./(1+(b01.*P.*exp(delU1./(8.314.*T)))) ...
                            + qs2.*(b02.*P.*exp(delU2./(8.314.*T)))./(1+(b02.*P.*exp(delU2./(8.314.*T))));
                    case 'DSS'
                        qvals(jj,kk) = qs1.*((b01.*P.*exp(delU1./(8.314.*T)))^gamma)./(1+(b01.*P.*exp(delU1./(8.314.*T)))^gamma) ...
                            + qs2.*((b02.*P.*exp(delU2./(8.314.*T)))^gamma)./(1+(b02.*P.*exp(delU2./(8.314.*T)))^gamma);
                    case 'SSS'
                        qvals(jj,kk) = qs1.*((b01.*P.*exp(delU1./(8.314.*T)))^gamma)./(1+(b01.*P.*exp(delU1./(8.314.*T)))^gamma) ...
                            + qs2.*((b02.*P.*exp(delU2./(8.314.*T)))^gamma)./(1+(b02.*P.*exp(delU2./(8.314.*T)))^gamma);
                end
            end
        end
        if ~flagConcUnits
            figure(1)
            if flagStats
                subplot(1,3,1)
            else
            end
            scatter(uncBounds(:,1),uncBounds(:,2),0.5,'MarkerEdgeColor','b','MarkerEdgeAlpha',0.2)
            hold on
            scatter(uncBounds(:,1),uncBounds(:,3),0.5,'MarkerEdgeColor','b','MarkerEdgeAlpha',0.2)
            for kk = 1:length(Tvals)
                plot(Pvals,qvals(:,kk),'-k','LineWidth',1.5);
            end
            outFit = [Pvals' qvals];
            plot(x,z,'ob');
            xlabel('Pressure [bar]');
            ylabel('Amount adsorbed [mol/kg]');
            % quantile-quantile plot of experimental data vs normal distribution
            box on
            set(gca,'YScale','linear','XScale','linear','FontSize',15,'LineWidth',1)
            grid on; axis square
            if flagStats
                subplot(1,3,2)
                pd = fitdist(z-qfit,'Normal');
                qqplot(qfit,pd);
                box on
                title('')
                set(gca,'YScale','linear','XScale','linear','FontSize',15,'LineWidth',1)
                grid on; axis square
                subplot(1,3,3)
                % Histogram for the error (experimental - fitted q) overlayed by the normal
                % distribution fitted for this error
                x_values = linspace(min(z-qfit),max(z-qfit),100);
                y_values = pdf(pd,x_values);
                yyaxis right
                plot(x_values,y_values,'LineWidth',0.5)
                ylabel('PDF [-]');
                yyaxis left
                histogram(z-qfit,15);
                xlabel('error [exp - model]');
                ylabel('Number of points, N_t [-]');
                box on
                set(gca,'YScale','linear','XScale','linear','FontSize',15,'LineWidth',1)
                grid on; axis square
            else
            end
            if ~flagContour
            else
                % Plot the contour plots for the objective function (MLE or WSS) around the
                % optimal parameter values in the combinations qs1-qs2, b01-delU1, and
                % b02-delU2.
                switch isothermModel
                    case 'DSL'
                        if length(unique(y)) == 1
                        else
                            Z = generateObjfunContour(x,y,z,nbins,isothermModel,fittingMethod, isoRef,parameters);
                        end
                    case 'DSS'
                        if length(unique(y)) == 1
                        else
                            Z = generateObjfunContour(x,y,z,nbins,isothermModel,fittingMethod, isoRef,parameters);
                        end
                end
            end
        end
end

% Save outputdata
if flagConcUnits
    isothermData.experiment = [x./(1e5./(8.314.*y)) z y];
else
    isothermData.experiment = [x z y];
end
headerRow = [NaN unique(y)'];
switch isothermModel
    case 'VIRIAL'
        isothermData.isothermFit = [headerRow;qvals(1,:)' lnPvals];
        isothermData.experiment = [z log(x) y];
    otherwise
        isothermData.isothermFit = [headerRow;Pvals(1,:)' qvals];
end
isothermData.confidenceRegion = outScatter;
isothermData.confidenceBounds = uncBounds;
isothermData.isothermParameters = [parameters' conRange95Disp];
isothermData.gitCommitID = gitCommitID;
if flagConcUnits
    isothermData.isothermFit(2:end,1) = linspace(0,x(find(x==max(max(x))))./(1e5./(8.314.*y(find(x==max(max(x)))))),200)';
    isothermData.confidenceBounds(:,1) = uncBounds(:,1)./(1e5./(8.314.*uncBounds(find(x==max(max(x))),4)));
end



if ~saveFlag
else
    filename = input('Enter file name: ','s');
    currentDate=datestr(date,'mmddyy');
    if exist(['..',filesep,'IsothermFittingTool',filesep','fittingResults'],'dir') == 7
        % Save the fitting results for further use
        save(['..',filesep,'IsothermFittingTool',filesep','fittingResults',filesep,filename,'_',currentDate],'isothermData');
    else
        % Create the fitting results folder if it does not exist
        mkdir(['..',filesep,'IsothermFittingTool',filesep','fittingResults'])
        % Save the calibration data for further use
        save(['..',filesep,'IsothermFittingTool',filesep','fittingResults',filesep,filename,'_',currentDate],'isothermData');
    end
end

if flagConcUnits
    figure
    plot(isothermData.isothermFit(:,1),isothermData.isothermFit(:,2:end),'-k','LineWidth',1.5)
    hold on;
    plot(isothermData.experiment(:,1),isothermData.experiment(:,2),'ok')
    scatter(isothermData.confidenceBounds(:,1),isothermData.confidenceBounds(:,2),0.5,'MarkerEdgeColor','b','MarkerEdgeAlpha',0.5);
    scatter(isothermData.confidenceBounds(:,1),isothermData.confidenceBounds(:,3),0.5,'MarkerEdgeColor','b','MarkerEdgeAlpha',0.5);
    xlabel('Concentration [mol/m3]');
    ylabel('Adsorbed amount [mol/kg]');
end


if strcmp(currentDir(end),'ERASE')
    cd ..
end