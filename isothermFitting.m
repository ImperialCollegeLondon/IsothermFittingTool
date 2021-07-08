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
clc; clear all;
% For new data create a new variable (HOME -> New Variable) named 'fitData'
% (rename on workspace) and save the data as *.mat file in 3 column format
% with Pressure (bar), adsorbed amount (-), temperature (K) respectively
% with any name of your choice.
% RUN THIS SCRIPT.
% Load input experimental data from *.mat or *.csv file via prompt
uiopen
% Determine number of bins you want the experimental data to be binned to
% in terms of the total pressure range
% (for Weighted sum of squares method ONLY).
% IF ERROR --> Reduce number of bins until error is gone
nbins = 1;
% Select isotherm model for fitting
% DSL = Dual site Langmuir. SSL = Single site Langmuir. DSS = Dual site
% Sips. SSS = Single site Sips
isothermModel = 'DSS';
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
flagFixQsat = 1;
% Flag for saving output in a matfile
saveFlag = 0;
% IF YOU ARE FIXING SATURATION CAPACITY ENTER THE CO2 SATURATION CAPACITIES
% FOR THE RELEVANT MODEL BELOW
qs1 = 3.6749e+00;
qs2 = 6.6845e+00;
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
% for concentration units, convert p to concentration
if ~flagFixQsat
    if flagConcUnits
        x = 1e5*x./(8.314.*y);
    end
    rng default % For reproducibility
    % Create gs, a GlobalSearch solver with its properties set to the defaults.
    gs = GlobalSearch('NumTrialPoints',1400,'NumStageOnePoints',1000,'Display','off');
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
                        optfunc = @(par) generateMLEfun(x, y, z, nbins, 'DSL', par(1), par(2), par(3), ...
                            par(4), 0, 0);
                    else
                        % Generate objective function for MLE method
                        optfunc = @(par) generateMLEfun(x, y, z, nbins, 'DSL', par(1), par(2), par(3), ...
                            par(4), par(5), par(6));
                    end
                case 'WSS'
                    % Generate objective function for WSS method
                    if length(unique(y)) == 1
                        % Generate objective function for MLE method
                        optfunc = @(par) generateWSSfun(x, y, z, nbins, 'DSL', par(1), par(2), par(3), ...
                            par(4), 0, 0);
                    else
                        % Generate objective function for MLE method
                        optfunc = @(par) generateWSSfun(x, y, z, nbins, 'DSL', par(1), par(2), par(3), ...
                            par(4), par(5), par(6));
                    end
            end
            % Initial conditions, lower bounds, and upper bounds for parameters
            % in DSL isotherm model
            if length(unique(y)) == 1 % if only one temperature
                x0 = [3, 3,1e-5,1e-5];
                lb = [0,0,0,0];
                ub = [20,20,10,10];
            else
                x0 = [3,3,1e-5,1e-5,5e4,5e4];
                lb = [0,0,0,0,0,0];
                ub = [20,20,1,1,7e4,7e4];
            end
            % Create global optimisation problem with solver 'fmincon' and
            % other bounds
            problem = createOptimProblem('fmincon','x0',x0,'objective',optfunc,'lb',lb,'ub',ub);
            % Solve the optimisation problem to obtain the isotherm parameters
            % for the fit
            [parVals, fval]= run(gs,problem);
            % Set fitted parameter values for isotherm model calculation
            qs1   = parVals(1);
            qs2   = parVals(2);
            b01   = parVals(3);
            b02   = parVals(4);
            if length(unique(y)) == 1 % if only one temperature
                delU1 = 0;
                delU2 = 0;
            else
                delU1 = parVals(5);
                delU2 = parVals(6);
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
                        optfunc = @(par) generateMLEfun(x, y, z, nbins, 'DSL', par(1), 0, par(2), ...
                            0, 0, 0);
                    else
                        % Generate objective function for MLE method
                        optfunc = @(par) generateMLEfun(x, y, z, nbins, 'DSL', par(1), 0, par(2), ...
                            0, par(3), 0);
                    end
                case 'WSS'
                    % Generate objective function for WSS method
                    if length(unique(y)) == 1
                        % Generate objective function for MLE method
                        optfunc = @(par) generateWSSfun(x, y, z, nbins, 'DSL', par(1), 0, par(2), ...
                            0, 0, 0);
                    else
                        % Generate objective function for MLE method
                        optfunc = @(par) generateWSSfun(x, y, z, nbins, 'DSL', par(1), 0, par(2), ...
                            0, par(3), 0);
                    end
            end
            % Initial conditions, lower bounds, and upper bounds for parameters
            % in DSL isotherm model
            if length(unique(y)) == 1 % if only one temperature
                x0 = [3,1e-5];
                lb = [0,0];
                ub = [20,10];
            else
                x0 = [3,1e-5,1e4];
                lb = [0,0,0];
                ub = [20,1,5e4];
            end
            % Create global optimisation problem with solver 'fmincon' and
            % other bounds
            problem = createOptimProblem('fmincon','x0',x0,'objective',optfunc,'lb',lb,'ub',ub);
            % Solve the optimisation problem to obtain the isotherm parameters
            % for the fit
            [parVals, fval]= run(gs,problem);
            % Set fitted parameter values for isotherm model calculation
            qs1   = parVals(1);
            qs2   = 0;
            b01   = parVals(2);
            b02   = 0;
            if length(unique(y)) == 1 % if only one temperature
                delU1 = 0;
            else
                delU1 = parVals(3);
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
                    optfunc = @(par) generateMLEfun(x, y, z, nbins, isothermModel, par(1), par(2), par(3), ...
                        par(4), par(5), par(6), par(7));
                case 'WSS'
                    % Generate objective function for WSS method
                    optfunc = @(par) generateWSSfun(x, y, z, nbins, isothermModel, par(1), par(2), par(3), ...
                        par(4), par(5), par(6), par(7));
            end
            % Initial conditions, lower bounds, and upper bounds for parameters
            % in DSL isotherm model
            x0 = [3,3,1e-5,1e-5,2e4,2e4,1];
            lb = [0,0,0,0,0,0,0];
            ub = [20,20,1,1,8e4,8e4,2];
            % Create global optimisation problem with solver 'fmincon' and
            % other bounds
            problem = createOptimProblem('fmincon','x0',x0,'objective',optfunc,'lb',lb,'ub',ub);
            % Solve the optimisation problem to obtain the isotherm parameters
            % for the fit
            [parVals, fval]= run(gs,problem);
            % Set fitted parameter values for isotherm model calculation
            qs1   = parVals(1);
            qs2   = parVals(2);
            b01   = parVals(3);
            b02   = parVals(4);
            delU1 = parVals(5);
            delU2 = parVals(6);
            gamma = parVals(7);
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
                    optfunc = @(par) generateMLEfun(x, y, z, nbins, 'DSS', par(1), 0, par(2), ...
                        0, par(3), 0, par(4));
                case 'WSS'
                    % Generate objective function for WSS method
                    optfunc = @(par) generateWSSfun(x, y, z, nbins, 'DSS', par(1), 0, par(2), ...
                        0, par(3), 0, par(4));
            end
            % Initial conditions, lower bounds, and upper bounds for parameters
            % in DSL isotherm model
            x0 = [3,1e-5,1e4,1];
            lb = [0,0,0,0];
            ub = [20,1,9e4,2];
            % Create global optimisation problem with solver 'fmincon' and
            % other bounds
            problem = createOptimProblem('fmincon','x0',x0,'objective',optfunc,'lb',lb,'ub',ub);
            % Solve the optimisation problem to obtain the isotherm parameters
            % for the fit
            [parVals, fval]= run(gs,problem);
            % Set fitted parameter values for isotherm model calculation
            qs1   = parVals(1);
            qs2   = 0;
            b01   = parVals(2);
            b02   = 0;
            delU1 = parVals(3);
            delU2 = 0;
            gamma = parVals(4);
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
                        optfunc = @(par) generateMLEfun(x, y, z, nbins, 'DSL', qs1, qs2, par(1), ...
                            par(2), 0, 0);
                    else
                        % Generate objective function for MLE method
                        optfunc = @(par) generateMLEfun(x, y, z, nbins, 'DSL', qs1, qs2, par(1), ...
                            par(2), par(3), par(4));
                    end
                case 'WSS'
                    % Generate objective function for WSS method
                    if length(unique(y)) == 1
                        % Generate objective function for MLE method
                        optfunc = @(par) generateWSSfun(x, y, z, nbins, 'DSL', qs1, qs2, par(1), ...
                            par(2), 0, 0);
                    else
                        % Generate objective function for MLE method
                        optfunc = @(par) generateWSSfun(x, y, z, nbins, 'DSL', qs1, qs2, par(1), ...
                            par(2), par(3), par(4));
                    end
            end
            % Initial conditions, lower bounds, and upper bounds for parameters
            % in DSL isotherm model
            if length(unique(y)) == 1 % if only one temperature
                x0 = [1e-5,1e-5];
                lb = [0,0,0,0];
                ub = [20,20,10,10];
            else
                x0 = [1e-5,1e-5,1e4,1e4];
                lb = [0,0,0,0];
                ub = [1,1,7e4,7e4];
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
            b01   = parVals(1);
            b02   = parVals(2);
            if length(unique(y)) == 1 % if only one temperature
                delU1 = 0;
                delU2 = 0;
            else
                delU1 = parVals(3);
                delU2 = parVals(4);
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
                        optfunc = @(par) generateMLEfun(x, y, z, nbins, 'DSL', qs1, 0, par(1), ...
                            0, 0, 0);
                    else
                        % Generate objective function for MLE method
                        optfunc = @(par) generateMLEfun(x, y, z, nbins, 'DSL', qs1, 0, par(1), ...
                            0, par(2), 0);
                    end
                case 'WSS'
                    % Generate objective function for WSS method
                    if length(unique(y)) == 1
                        % Generate objective function for MLE method
                        optfunc = @(par) generateWSSfun(x, y, z, nbins, 'DSL', qs1, 0, par(1), ...
                            0, 0, 0);
                    else
                        % Generate objective function for MLE method
                        optfunc = @(par) generateWSSfun(x, y, z, nbins, 'DSL', qs1, 0, par(1), ...
                            0, par(2), 0);
                    end
            end
            % Initial conditions, lower bounds, and upper bounds for parameters
            % in DSL isotherm model
            if length(unique(y)) == 1 % if only one temperature
                x0 = [1e-5];
                lb = [0];
                ub = [10];
            else
                x0 = [1e-5,3e4];
                lb = [0,0];
                ub = [1,5e4];
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
            b01   = parVals(1);
            b02   = 0;
            if length(unique(y)) == 1 % if only one temperature
                delU1 = 0;
            else
                delU1 = parVals(2);
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
                    optfunc = @(par) generateMLEfun(x, y, z, nbins, isothermModel, qs1, qs2, par(1), ...
                        par(2), par(3), par(4), par(5));
                case 'WSS'
                    % Generate objective function for WSS method
                    optfunc = @(par) generateWSSfun(x, y, z, nbins, isothermModel, qs1, qs2, par(1), ...
                        par(2), par(3), par(4), par(5));
            end
            % Initial conditions, lower bounds, and upper bounds for parameters
            % in DSL isotherm model
            x0 = [1e-5,1e-5,3e4,3e4,1];
            lb = [0,0,0,0,0];
            ub = [1,1,5e4,5e4,2];
            % Create global optimisation problem with solver 'fmincon' and
            % other bounds
            problem = createOptimProblem('fmincon','x0',x0,'objective',optfunc,'lb',lb,'ub',ub);
            % Solve the optimisation problem to obtain the isotherm parameters
            % for the fit
            [parVals, fval]= run(gs,problem);
            % Set fitted parameter values for isotherm model calculation
            qs1   = qs1;
            qs2   = qs2;
            b01   = parVals(1);
            b02   = parVals(2);
            delU1 = parVals(3);
            delU2 = parVals(4);
            gamma = parVals(5);
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
                    optfunc = @(par) generateMLEfun(x, y, z, nbins, 'DSS', qs1, 0, par(1), ...
                        0, par(2), 0, par(3));
                case 'WSS'
                    % Generate objective function for WSS method
                    optfunc = @(par) generateWSSfun(x, y, z, nbins, 'DSS', qs1, 0, par(1), ...
                        0, par(2), 0, par(3));
            end
            % Initial conditions, lower bounds, and upper bounds for parameters
            % in DSL isotherm model
            x0 = [1e-5,3e4,1];
            lb = [0,0,0];
            ub = [1,5e4,2];
            % Create global optimisation problem with solver 'fmincon' and
            % other bounds
            problem = createOptimProblem('fmincon','x0',x0,'objective',optfunc,'lb',lb,'ub',ub);
            % Solve the optimisation problem to obtain the isotherm parameters
            % for the fit
            [parVals, fval]= run(gs,problem);
            % Set fitted parameter values for isotherm model calculation
            qs1   = qs1;
            qs2   = 0;
            b01   = parVals(1);
            b02   = 0;
            delU1 = parVals(2);
            delU2 = 0;
            gamma = parVals(3);
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
        [outScatter]=generateUncertaintySpread(x,y,'DSL',parameters,conRange95);
    case 'SSL'
        [outScatter]=generateUncertaintySpread(x,y,'DSL',parameters,conRange95);
    case 'DSS'
        [outScatter]=generateUncertaintySpread(x,y,'DSS',parameters,conRange95);
    case 'SSS'
        [outScatter]=generateUncertaintySpread(x,y,'DSS',parameters,conRange95);
end
% plot of experimental data and fitted data (q vs P)
Pvals = linspace(0,max(x),200);
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
if flagConcUnits
    Pvals = Pvals./(1e5./(8.314.*y));
    x = x./(1e5./(8.314.*y));
    outScatter(:,1) = outScatter(:,1)./(1e5./(8.314.*outScatter(:,3)));
end
if ~flagConcUnits
%     figure
    if flagStats
        subplot(1,3,1)
    else
    end
    scatter(outScatter(:,1),outScatter(:,2),5,'MarkerFaceColor','#D9E9FC','MarkerEdgeAlpha',0.05)
    hold on
    for kk = 1:length(Tvals)
        plot(Pvals,qvals(:,kk),'-k','LineWidth',1.5);
    end
    outFit = [Pvals' qvals];
    plot(x,z,'ok');
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
                    Z = generateObjfunContour(x,y,z,nbins,isothermModel,fittingMethod,parameters);
                end
            case 'DSS'
                if length(unique(y)) == 1
                else
                    Z = generateObjfunContour(x,y,z,nbins,isothermModel,fittingMethod,parameters);
                end
        end
    end
end

% Save outputdata
isothermData.experiment = fitData;
isothermData.isothermFit = [x qfit y];
isothermData.confidenceRegion = outScatter;
isothermData.isothermParameters = [parameters' conRange95Disp];
if ~saveFlag
else
    filename = input('Enter file name: ','s');
    currentDate=datestr(date,'mmddyy');
    if exist(['..',filesep,'IsothermFittingTool',filesep','fittingResults'],'dir') == 7
        % Save the fitting results for further use
        save(['..',filesep,'IsothermFittingTool',filesep','fittingResults',filesep,filename,'_',currentDate],'isothermData');
        clear all
    else
        % Create the fitting results folder if it does not exist
        mkdir(['..',filesep,'IsothermFittingTool',filesep','fittingResults'])
        % Save the calibration data for further use
        save(['..',filesep,'IsothermFittingTool',filesep','fittingResults',filesep,filename,'_',currentDate],'isothermData');
        clear all
    end
end
