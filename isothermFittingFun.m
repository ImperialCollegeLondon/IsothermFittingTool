function [varargout] = isothermFittingFun(isothermModel, flagConcUnits, saveFlag, fitData)
%% GENERATING AND SOLVING GLOBAL OPTIMISATION PROBLEM
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
% Determine number of bins you want the experimental data to be binned to
% in terms of the total pressure range
% (for Weighted sum of squares method ONLY).
% IF ERROR --> Reduce number of bins until error is gone
nbins = 1;
% Select fitting method.
% WSS = weighted sum of squares, MLE = maximum log-likelihood estimator
% MLE is preferred for data that is from a single source where the error is
% likely to be normally distributed with a mean of 0.
% WSS is preferred for fitting data for cases where the error might not be
% random and not be normally distributed (data from different sources)
fittingMethod = 'MLE';
% Flag for plotting objective function contour plots for dual site models
flagContour = 0;
% Flag for plotting statistical plots (q-q plot and error distribution)
flagStats = 0;
% Flag for fixing saturation capacities (0 for CO2 fitting, 1 for other
% gases)
flagFixQsat = 0;
% IF YOU ARE FIXING SATURATION CAPACITY ENTER THE CO2 SATURATION CAPACITIES
% FOR THE RELEVANT MODEL BELOW
qs1 = 1.9834e+00;
qs2 = 8.9856e+00;
% Pressure (x), adsorbed amount (z), Temperature (y) data from input
x = fitData(:,1);
z = fitData(:,2);
y = fitData(:,3);
% Reference isotherm parameters for non-dimensionalisation [qs1 qs2 b01 b02 delU1 delU2]
refValsP = [10,10,1e-1,1e-1,5e4,5e4];
refValsC = [10,10,1e-5,1e-5,5e4,5e4];
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
    case  'TSL'
        % Reference isotherm parameters for non-dimensionalisation
        if flagConcUnits
            isoRef = [refValsC refValsC([1 3 5])];
        else
            isoRef = [refValsP refValsP([1 3 5])];
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
    case  'TOTH'
        % Reference isotherm parameters for non-dimensionalisation
        if flagConcUnits
            isoRef = [refValsC 2];
        else
            isoRef = [refValsP 2];
        end
    case  'HDSL'
        % Reference isotherm parameters for non-dimensionalisation
        if flagConcUnits
            isoRef = [refValsC refValsC([1 3 5])];
        else
            isoRef = [refValsP refValsP([1 3 5])];
        end
    case  'HSSL'
        % Reference isotherm parameters for non-dimensionalisation
        if flagConcUnits
            isoRef = [refValsC refValsC([1 3 5])];
        else
            isoRef = [refValsP refValsP([1 3 5])];
        end
end
% for concentration units, convert pressure to concentration
if ~flagFixQsat
    if flagConcUnits
        x = 1e5*x./(8.314.*y);
    end
    rng default % For reproducibility
    % Create gs, a GlobalSearch solver with its properties set to the defaults.
    gs = GlobalSearch('NumTrialPoints',2000,'NumStageOnePoints',200,'Display','off');
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
                x0 = [0.5,0.5,0.5,0.5];
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
            [conRange95] = conrangeEllipse(x, y, z, qfit,fittingMethod,isoRef, 'DSL', qs1, qs2, b01, b02, delU1, delU2);
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
                        fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
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
                        fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
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
            [conRange95] = conrangeEllipse(x, y, z, qfit,fittingMethod,isoRef, 'DSL', qs1, qs2, b01, b02, delU1, delU2);
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
                        fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
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
                        fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            end
        case 'TSL'
            % Set objective function based on fitting method
            switch fittingMethod
                case 'MLE'
                    % Number of bins is automatically set to 1 for MLE as
                    % weights cannot be assigned in MLE
                    nbins =1;
                    if length(unique(y)) == 1 % if only one temperature
                        % Generate objective function for MLE method
                        optfunc = @(par) generateMLEfun(x, y, z, nbins, 'TSL', isoRef, par(1), par(2), par(3), ...
                            par(4), 0, 0, par(5), par(6), 0);
                    else
                        % Generate objective function for MLE method
                        optfunc = @(par) generateMLEfun(x, y, z, nbins, 'TSL', isoRef, par(1), par(2), par(3), ...
                            par(4), par(5), par(6), par(7), par(8),  par(9));
                    end
                case 'WSS'
                    % Generate objective function for WSS method
                    if length(unique(y)) == 1
                        % Generate objective function for MLE method
                        optfunc = @(par) generateWSSfun(x, y, z, nbins, 'TSL', isoRef, par(1), par(2), par(3), ...
                            par(4), 0, 0, par(5), par(6),  0);
                    else
                        % Generate objective function for MLE method
                        optfunc = @(par) generateWSSfun(x, y, z, nbins, 'TSL', isoRef, par(1), par(2), par(3), ...
                            par(4), par(5), par(6), par(7), par(8),  par(9));
                    end
            end
            % Initial conditions, lower bounds, and upper bounds for parameters
            % in DSL isotherm model
            if length(unique(y)) == 1 % if only one temperature
                x0 = [0.5,0.5,0.5,0.5,0.5,0.5];
                lb = [0,0,0,0,0,0];
                ub = [1,1,1,1,1,1];
            else
                x0 = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5];
                lb = [0,0,0,0,0,0,0,0,0];
                ub = [1,1,1,1,1,1,1,1,1];
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
                qs3   = parVals(5).*isoRef(7);
                b03   = parVals(6).*isoRef(8);
                delU3 = 0;
            else
                delU1 = parVals(5).*isoRef(5);
                delU2 = parVals(6).*isoRef(6);
                qs3   = parVals(7).*isoRef(7);
                b03   = parVals(8).*isoRef(8);
                delU3 = parVals(9).*isoRef(9);
            end
            
            % Calculate fitted isotherm loadings for conditions (P,T)
            % corresponding to experimental data
            qfit  = qs1.*(b01.*x.*exp(delU1./(8.314.*y)))./(1+(b01.*x.*exp(delU1./(8.314.*y)))) ...
                + qs2.*(b02.*x.*exp(delU2./(8.314.*y)))./(1+(b02.*x.*exp(delU2./(8.314.*y)))) ...
                + qs3.*(b03.*x.*exp(delU3./(8.314.*y)))./(1+(b03.*x.*exp(delU3./(8.314.*y))));
            % Calculate ellipsoidal confidence intervals (delta parameter) for
            % fitted parameters
            parameters = [qs1, qs2, b01, b02, delU1, delU2, qs3, b03, delU3];
            parameters(isnan(parameters))=0;
            [conRange95] = conrangeEllipse(x, y, z, qfit,fittingMethod,isoRef, 'TSL', qs1, qs2, b01, b02, delU1, delU2, qs3, b03, delU3);
            conRange95(isnan(conRange95))=0;
            % Convert confidence intervals to percentage error
            percentageError = conRange95./parameters' *100;
            fprintf('Isotherm model: %s \n', isothermModel);
            if ~flagConcUnits
                parNames = ["qs1" "qs2" "qs3" "b01" "b02" "b03" "delU1" "delU2" "delU3"];
                units = ["mol/kg" "mol/kg" "mol/kg" "1/bar" "1/bar" "1/bar" "J/mol" "J/mol" "J/mol"];
                parsDisp = [qs1 qs2 qs3 b01 b02 b03 delU1 delU2 delU3];
                conRange95Disp = [conRange95(1) conRange95(2) conRange95(7) conRange95(3) conRange95(4) conRange95(8) conRange95(5) conRange95(6) conRange95(9)]';
                for ii = 1:length(parameters)
                    if parsDisp(ii) == 0
                    else
                        fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            else
                parNames = ["qs1" "qs2" "qs3" "b01" "b02" "b03" "delU1" "delU2" "delU3"];
                units = ["mol/kg" "mol/kg" "mol/kg" "1/bar" "1/bar" "1/bar" "J/mol" "J/mol" "J/mol"];
                parsDisp = [qs1 qs2 qs3 b01 b02 b03 delU1 delU2 delU3];
                conRange95Disp = [conRange95(1) conRange95(2) conRange95(7) conRange95(3) conRange95(4) conRange95(8) conRange95(5) conRange95(6) conRange95(9)]';
                for ii = 1:length(parsDisp)
                    if parsDisp(ii) == 0
                    else
                        fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            end
        case 'HDSL'
            % Set objective function based on fitting method
            switch fittingMethod
                case 'MLE'
                    % Number of bins is automatically set to 1 for MLE as
                    % weights cannot be assigned in MLE
                    nbins =1;
                    if length(unique(y)) == 1 % if only one temperature
                        % Generate objective function for MLE method
                        optfunc = @(par) generateMLEfun(x, y, z, nbins, 'HDSL', isoRef, par(1), par(2), par(3), ...
                            par(4), 0, 0, par(5), par(6), 0);
                    else
                        % Generate objective function for MLE method
                        optfunc = @(par) generateMLEfun(x, y, z, nbins, 'HDSL', isoRef, par(1), par(2), par(3), ...
                            par(4), par(5), par(6), par(7), par(8),  par(9));
                    end
                case 'WSS'
                    % Generate objective function for WSS method
                    if length(unique(y)) == 1
                        % Generate objective function for MLE method
                        optfunc = @(par) generateWSSfun(x, y, z, nbins, 'HDSL', isoRef, par(1), par(2), par(3), ...
                            par(4), 0, 0, par(5), par(6),  0);
                    else
                        % Generate objective function for MLE method
                        optfunc = @(par) generateWSSfun(x, y, z, nbins, 'HDSL', isoRef, par(1), par(2), par(3), ...
                            par(4), par(5), par(6), par(7), par(8),  par(9));
                    end
            end
            % Initial conditions, lower bounds, and upper bounds for parameters
            % in DSL isotherm model
            if length(unique(y)) == 1 % if only one temperature
                x0 = [0.5,0.5,0.5,0.5,0.5,0.5];
                lb = [0,0,0,0,0,0];
                ub = [1,1,1,1,1,1];
            else
                x0 = [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5];
                lb = [0,0,0,0,0,0,0,0,0];
                ub = [1,1,1,1,1,1,1,1,1];
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
                qsH   = parVals(5).*isoRef(7);
                b0H   = parVals(6).*isoRef(8);
                delUH = 0;
            else
                delU1 = parVals(5).*isoRef(5);
                delU2 = parVals(6).*isoRef(6);
                qsH   = parVals(7).*isoRef(7);
                b0H   = parVals(8).*isoRef(8);
                delUH = parVals(9).*isoRef(9);
            end
            
            % Calculate fitted isotherm loadings for conditions (P,T)
            % corresponding to experimental data
            qfit  = qs1.*(b01.*x.*exp(delU1./(8.314.*y)))./(1+(b01.*x.*exp(delU1./(8.314.*y)))) ...
                + qs2.*(b02.*x.*exp(delU2./(8.314.*y)))./(1+(b02.*x.*exp(delU2./(8.314.*y)))) ...
                + qsH.*(b0H.*x.*exp(delUH./(8.314.*y)));
            % Calculate ellipsoidal confidence intervals (delta parameter) for
            % fitted parameters
            parameters = [qs1, qs2, b01, b02, delU1, delU2, qsH, b0H, delUH];
            parameters(isnan(parameters))=0;
            [conRange95] = conrangeEllipse(x, y, z, qfit,fittingMethod,isoRef, 'HDSL', qs1, qs2, b01, b02, delU1, delU2, qsH, b0H, delUH);
            conRange95(isnan(conRange95))=0;
            % Convert confidence intervals to percentage error
            percentageError = conRange95./parameters' *100;
            fprintf('Isotherm model: %s \n', isothermModel);
            if ~flagConcUnits
                parNames = ["qs1" "qs2" "b01" "b02" "delU1" "delU2" "qsH" "b0H" "delUH" ];
                units = ["mol/kg" "mol/kg" "1/bar" "1/bar" "J/mol" "J/mol" "mol/kg" "1/bar" "J/mol"];
                parsDisp = parameters;
                conRange95Disp = conRange95;
                for ii = 1:length(parameters)
                    if parsDisp(ii) == 0
                    else
                        fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            else
                parNames = ["qs1" "qs2" "b01" "b02" "delU1" "delU2" "qsH" "b0H" "delUH" ];
                units = ["mol/kg" "mol/kg" "1/bar" "1/bar" "J/mol" "J/mol" "mol/kg" "1/bar" "J/mol"];
                parsDisp = parameters;
                conRange95Disp = conRange95;
                for ii = 1:length(parsDisp)
                    if parsDisp(ii) == 0
                    else
                        fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            end
        case 'HSSL'
            % Set objective function based on fitting method
            switch fittingMethod
                case 'MLE'
                    % Number of bins is automatically set to 1 for MLE as
                    % weights cannot be assigned in MLE
                    nbins =1;
                    if length(unique(y)) == 1 % if only one temperature
                        % Generate objective function for MLE method
                        optfunc = @(par) generateMLEfun(x, y, z, nbins, 'HDSL', isoRef, par(1), 0, par(2), ...
                            0, 0, 0, par(3), par(4), 0);
                    else
                        % Generate objective function for MLE method
                        optfunc = @(par) generateMLEfun(x, y, z, nbins, 'HDSL', isoRef, par(1), 0, par(2), ...
                            0, par(3), 0, par(4), par(5),  par(6));
                    end
                case 'WSS'
                    % Generate objective function for WSS method
                    if length(unique(y)) == 1 % if only one temperature
                        % Generate objective function for MLE method
                        optfunc = @(par) generateWSSfun(x, y, z, nbins, 'HDSL', isoRef, par(1), 0, par(2), ...
                            0, 0, 0, par(3), par(4), 0);
                    else
                        % Generate objective function for MLE method
                        optfunc = @(par) generateWSSfun(x, y, z, nbins, 'HDSL', isoRef, par(1), 0, par(2), ...
                            0, par(3), 0, par(4), par(5),  par(6));
                    end
            end
            % Initial conditions, lower bounds, and upper bounds for parameters
            % in DSL isotherm model
            if length(unique(y)) == 1 % if only one temperature
                x0 = [0.5,0.5,0.5,0.5];
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
            qs2   = 0;
            b01   = parVals(2).*isoRef(3);
            b02   = 0;
            if length(unique(y)) == 1 % if only one temperature
                delU1 = 0;
                delU2 = 0;
                qsH   = parVals(3).*isoRef(7);
                b0H   = parVals(4).*isoRef(8);
                delUH = 0;
            else
                delU1 = parVals(3).*isoRef(5);
                delU2 = 0;
                qsH   = parVals(4).*isoRef(7);
                b0H   = parVals(5).*isoRef(8);
                delUH = parVals(6).*isoRef(9);
            end
            
            % Calculate fitted isotherm loadings for conditions (P,T)
            % corresponding to experimental data
            qfit  = qs1.*(b01.*x.*exp(delU1./(8.314.*y)))./(1+(b01.*x.*exp(delU1./(8.314.*y)))) ...
                + qs2.*(b02.*x.*exp(delU2./(8.314.*y)))./(1+(b02.*x.*exp(delU2./(8.314.*y)))) ...
                + qsH.*(b0H.*x.*exp(delUH./(8.314.*y)));
            % Calculate ellipsoidal confidence intervals (delta parameter) for
            % fitted parameters
            parameters = [qs1, qs2, b01, b02, delU1, delU2, qsH, b0H, delUH];
            parameters(isnan(parameters))=0;
            [conRange95] = conrangeEllipse(x, y, z, qfit,fittingMethod,isoRef, 'HDSL', qs1, qs2, b01, b02, delU1, delU2, qsH, b0H, delUH);
            conRange95(isnan(conRange95))=0;
            % Convert confidence intervals to percentage error
            percentageError = conRange95./parameters' *100;
            fprintf('Isotherm model: %s \n', isothermModel);
            if ~flagConcUnits
                parNames = ["qs1" "qs2" "b01" "b02" "delU1" "delU2" "qsH" "b0H" "delUH" ];
                units = ["mol/kg" "mol/kg" "1/bar" "1/bar" "J/mol" "J/mol" "mol/kg" "1/bar" "J/mol"];
                parsDisp = parameters;
                conRange95Disp = conRange95;
                for ii = 1:length(parameters)
                    if parsDisp(ii) == 0
                    else
                        fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            else
                parNames = ["qs1" "qs2" "b01" "b02" "delU1" "delU2" "qsH" "b0H" "delUH" ];
                units = ["mol/kg" "mol/kg" "1/bar" "1/bar" "J/mol" "J/mol" "mol/kg" "1/bar" "J/mol"];
                parsDisp = parameters;
                conRange95Disp = conRange95;
                for ii = 1:length(parsDisp)
                    if parsDisp(ii) == 0
                    else
                        fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
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
            [conRange95] = conrangeEllipse(x, y, z, qfit,fittingMethod,isoRef, isothermModel, qs1, qs2, b01, b02, delU1, delU2, gamma);
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
                        fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
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
                        fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
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
            [conRange95] = conrangeEllipse(x, y, z, qfit,fittingMethod,isoRef, 'DSS', qs1, qs2, b01, b02, delU1, delU2, gamma);
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
                        fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
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
                        fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            end
        case 'TOTH'
            % Set objective function based on fitting method
            switch fittingMethod
                case 'MLE'
                    % Number of bins is automatically set to 1 for MLE as
                    % weights cannot be assigned in MLE
                    nbins =1;
                    % Generate objective function for MLE method
                    optfunc = @(par) generateMLEfun(x, y, z, nbins, 'TOTH', isoRef, par(1), 0, par(2), ...
                        0, par(3), 0, par(4));
                case 'WSS'
                    % Generate objective function for WSS method
                    optfunc = @(par) generateWSSfun(x, y, z, nbins, 'TOTH', isoRef, par(1), 0, par(2), ...
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
            toth = parVals(4).*isoRef(7);
            % Calculate fitted isotherm loadings for conditions (P,T)
            % corresponding to experimental data
            qfit  = qs1.*b01.*x.*exp(delU1./(8.314.*y))./(1+(b01.*x.*exp(delU1./(8.314.*y))).^toth).^(1./toth);
            % Calculate ellipsoidal confidence intervals (delta parameter) for
            % fitted parameters
            parameters = [qs1, qs2, b01, b02, delU1, delU2, toth];
            parameters(isnan(parameters))=0;
            [conRange95] = conrangeEllipse(x, y, z, qfit,fittingMethod,isoRef, 'TOTH', qs1, qs2, b01, b02, delU1, delU2, toth);
            conRange95(isnan(conRange95))=0;
            conRange95 = [conRange95(1) 0 conRange95(3) 0 conRange95(5) 0 conRange95(7)]';
            % Convert confidence intervals to percentage error
            percentageError = conRange95./parameters' *100;
            fprintf('Isotherm model: %s \n', isothermModel);
            if ~flagConcUnits
                parNames = ["qs1" "qs2" "b01" "b02" "delU1" "delU2" "toth"];
                units = ["mol/kg" "mol/kg" "1/bar" "1/bar" "J/mol" "J/mol" " "];
                parsDisp = parameters;
                conRange95Disp = conRange95;
                for ii = 1:length(parameters)
                    if parsDisp(ii) == 0
                    else
                        fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            else
                parNames = ["qs1" "qs2" "b01" "b02" "delU1" "delU2"  "toth"];
                units = ["mol/kg" "mol/kg" "m3/mol" "m3/mol" "J/mol" "J/mol"  " "];
                parsDisp = parameters;
                conRange95Disp = conRange95;
                for ii = 1:length(parsDisp)
                    if parsDisp(ii) == 0
                    else
                        fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            end
        case 'VIRIAL'
            % Reference isotherm parameters for non-dimensionalisation [qs1 qs2 b01 b02 delU1 delU2]
            refValsP = [10e3,10e3,10e3,10e3,10e3,10e3,10e3,10e3];
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
            lnPfit = [];
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
            [conRange95] = conrangeEllipse(x, y, z, lnPfit, fittingMethod,isoRef, 'VIRIAL', a0, a1, a2, a3, b0, b1, b2, b3);
            conRange95(isnan(conRange95))=0;
            conRange95 = [conRange95(1) conRange95(2) conRange95(3) conRange95(4) conRange95(5) conRange95(6) 0 0]';
            % Convert confidence intervals to percentage error
            percentageError = conRange95./parameters' *100;
            fprintf('Isotherm model: %s \n', isothermModel);
            parNames = ["a0" "a1" "a2" "a3" "b0" "b1" "b2" "b3"];
            units = ["K" "K/mol" "K/mol^2" "K/mol^3" " " "1/mol" "1/mol^2" "1/mol^3"];
            parsDisp = parameters;
            conRange95Disp = conRange95;
            for ii = 1:length(parameters)
                if parsDisp(ii) == 0
                else
                    fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                end
            end
            
        case 'VIRIAL2'
            % Reference isotherm parameters for non-dimensionalisation [qs1 qs2 b01 b02 delU1 delU2]
            refValsP = [10e3,10e3,10e3,10e3,10e3,10e3,10e3,10e3];
            refValsC = refValsP;
            isoRef = refValsC;
            % Set objective function based on fitting method
            switch fittingMethod
                case 'MLE'
                    % Number of bins is automatically set to 1 for MLE as
                    % weights cannot be assigned in MLE
                    nbins =1;
                    % Generate objective function for MLE method
                    optfunc = @(par) generateMLEfun(x, y, z, nbins, 'VIRIAL2', isoRef, par(1), par(2), ...
                        par(3), par(4), par(5), par(6), par(7), par(8));
                case 'WSS'
                    % Generate objective function for WSS method
                    optfunc = @(par) generateWSSfun(x, y, z, nbins, 'VIRIAL2', isoRef, par(1), par(2), ...
                        par(3), par(4), par(5), par(6), par(7), par(8));
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
            b2   = parVals(7).*isoRef(7);
            b3   = parVals(8).*isoRef(8);
            parameters = [a0, a1, a2, a3, b0, b1, b2, b3];
            parameters(isnan(parameters))=0;
            % Calculate fitted isotherm loadings for conditions (P,T)
            % corresponding to experimental data
            lnPfit = [];
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
            [conRange95] = conrangeEllipse(x, y, z, lnPfit, fittingMethod,isoRef, 'VIRIAL2', a0, a1, a2, a3, b0, b1, b2, b3);
            conRange95(isnan(conRange95))=0;
            conRange95 = [conRange95(1) conRange95(2) conRange95(3) conRange95(4) conRange95(5) conRange95(6) conRange95(7) conRange95(8)]';
            % Convert confidence intervals to percentage error
            percentageError = conRange95./parameters' *100;
            fprintf('Isotherm model: %s \n', isothermModel);
            parNames = ["a0" "a1" "a2" "a3" "b0" "b1" "b2" "b3"];
            units = ["K" "K/mol" "K/mol^2" "K/mol^3" " " "1/mol" "1/mol^2" "1/mol^3"];
            parsDisp = parameters;
            conRange95Disp = conRange95;
            for ii = 1:length(parameters)
                if parsDisp(ii) == 0
                else
                    fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
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
            [conRange95] = conrangeEllipse(x, y, z, qfit,fittingMethod,isoRef, 'DSL', qs1, qs2, b01, b02, delU1, delU2);
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
                        fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
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
                        fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
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
            [conRange95] = conrangeEllipse(x, y, z, qfit,fittingMethod,isoRef, 'DSL', qs1, qs2, b01, b02, delU1, delU2);
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
                        fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
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
                        fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
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
            [conRange95] = conrangeEllipse(x, y, z, qfit,fittingMethod,isoRef, isothermModel, qs1, qs2, b01, b02, delU1, delU2, gamma);
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
                        fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
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
                        fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
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
            [conRange95] = conrangeEllipse(x, y, z, qfit,fittingMethod,isoRef, 'DSS', qs1, qs2, b01, b02, delU1, delU2, gamma);
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
                        fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
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
                        fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            end
        case 'TOTH'
            % Set objective function based on fitting method
            switch fittingMethod
                case 'MLE'
                    % Number of bins is automatically set to 1 for MLE as
                    % weights cannot be assigned in MLE
                    nbins =1;
                    % Generate objective function for MLE method
                    optfunc = @(par) generateMLEfun(x, y, z, nbins, 'TOTH', isoRef, qs1, 0, par(1), ...
                        0, par(2), 0, par(3));
                case 'WSS'
                    % Generate objective function for WSS method
                    optfunc = @(par) generateWSSfun(x, y, z, nbins, 'TOTH', isoRef, qs1, 0, par(1), ...
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
            toth = parVals(3).*isoRef(7);
            % Calculate fitted isotherm loadings for conditions (P,T)
            % corresponding to experimental data
            qfit  = qs1.*b01.*x.*exp(delU1./(8.314.*y))./(1+(b01.*x.*exp(delU1./(8.314.*y))).^toth).^(1./toth);
            % Calculate ellipsoidal confidence intervals (delta parameter) for
            % fitted parameters
            parameters = [qs1, qs2, b01, b02, delU1, delU2, toth];
            parameters(isnan(parameters))=0;
            [conRange95] = conrangeEllipse(x, y, z, qfit,fittingMethod,isoRef, 'TOTH', qs1, qs2, b01, b02, delU1, delU2, toth);
            conRange95(isnan(conRange95))=0;
            % Convert confidence intervals to percentage error
            percentageError = conRange95./parameters' *100;
            fprintf('Isotherm model: %s \n', isothermModel);
            if ~flagConcUnits
                parNames = ["qs1" "qs2" "b01" "b02" "delU1" "delU2" "toth"];
                units = ["mol/kg" "mol/kg" "1/bar" "1/bar" "J/mol" "J/mol" " "];
                parsDisp = parameters;
                conRange95Disp = conRange95;
                for ii = 1:length(parameters)
                    if parsDisp(ii) == 0
                    else
                        fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
                    end
                end
            else
                parNames = ["qs1" "qs2" "b01" "b02" "delU1" "delU2"  "toth"];
                units = ["mol/kg" "mol/kg" "m3/mol" "m3/mol" "J/mol" "J/mol"  " "];
                parsDisp = parameters;
                conRange95Disp = conRange95;
                for ii = 1:length(parsDisp)
                    if parsDisp(ii) == 0
                    else
                        fprintf('%s = %5.4e ?? %5.4e %s \n',parNames(ii),parsDisp(ii),conRange95Disp(ii),units(ii));
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
    case 'TSL'
        [outScatter,uncBounds]=generateUncertaintySpread(x,y,z,'TSL',parameters,conRange95);
    case 'HDSL'
        [outScatter,uncBounds]=generateUncertaintySpread(x,y,z,'HDSL',parameters,conRange95);
    case 'HSSL'
        [outScatter,uncBounds]=generateUncertaintySpread(x,y,z,'HDSL',parameters,conRange95);
    case 'DSS'
        [outScatter,uncBounds]=generateUncertaintySpread(x,y,z,'DSS',parameters,conRange95);
    case 'SSS'
        [outScatter,uncBounds]=generateUncertaintySpread(x,y,z,'DSS',parameters,conRange95);
    case 'TOTH'
        [outScatter,uncBounds]=generateUncertaintySpread(x,y,z,'TOTH',parameters,conRange95);
    case 'VIRIAL'
        [outScatter,uncBounds]=generateUncertaintySpread(x,y,z,'VIRIAL',parameters,conRange95);
    case 'VIRIAL2'
        [outScatter,uncBounds]=generateUncertaintySpread(x,y,z,'VIRIAL2',parameters,conRange95);
end
switch isothermModel
    case {'VIRIAL','VIRIAL2'}
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
        xlim([0 max(x)*1.1]);
        ylim([0 max(z)*1.1]);
        box on
        set(gca,'YScale','linear','XScale','linear','FontSize',15,'LineWidth',1)
        grid on; axis square
        subplot(1,3,3)
        for kk = 1:length(Tvals)
            delH = -8.314.*(a0 +a1.*qvals+a2.*qvals.^2+a3.*qvals.^3)./1000;
            plot(qvals,delH,'-k','LineWidth',1.5);
        end
        xlabel('Amount adsorbed [mol/kg]');
        ylabel('-\Delta H_{ads} [kJ/mol]');
        ylim([0 max(delH)+5]);
        box on
        set(gca,'YScale','linear','XScale','linear','FontSize',15,'LineWidth',1)
        grid on; axis square
        set(gcf,'units','inch','position',[0,0,10,4])
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
                    case 'HDSL'
                        qvals(jj,kk) = qs1.*(b01.*P.*exp(delU1./(8.314.*T)))./(1+(b01.*P.*exp(delU1./(8.314.*T)))) ...
                            + qs2.*(b02.*P.*exp(delU2./(8.314.*T)))./(1+(b02.*P.*exp(delU2./(8.314.*T)))) ...
                            + qsH.*(b0H.*P.*exp(delUH./(8.314.*T)));
                    case 'TSL'
                        qvals(jj,kk) = qs1.*(b01.*P.*exp(delU1./(8.314.*T)))./(1+(b01.*P.*exp(delU1./(8.314.*T)))) ...
                            + qs2.*(b02.*P.*exp(delU2./(8.314.*T)))./(1+(b02.*P.*exp(delU2./(8.314.*T)))) ...
                            + qs3.*(b03.*P.*exp(delU3./(8.314.*T)))./(1+(b03.*P.*exp(delU3./(8.314.*T))));
                    case 'HSSL'
                        qvals(jj,kk) = qs1.*(b01.*P.*exp(delU1./(8.314.*T)))./(1+(b01.*P.*exp(delU1./(8.314.*T)))) ...
                            + qs2.*(b02.*P.*exp(delU2./(8.314.*T)))./(1+(b02.*P.*exp(delU2./(8.314.*T)))) ...
                            + qsH.*(b0H.*P.*exp(delUH./(8.314.*T)));
                    case 'DSS'
                        qvals(jj,kk) = qs1.*((b01.*P.*exp(delU1./(8.314.*T)))^gamma)./(1+(b01.*P.*exp(delU1./(8.314.*T)))^gamma) ...
                            + qs2.*((b02.*P.*exp(delU2./(8.314.*T)))^gamma)./(1+(b02.*P.*exp(delU2./(8.314.*T)))^gamma);
                    case 'SSS'
                        qvals(jj,kk) = qs1.*((b01.*P.*exp(delU1./(8.314.*T)))^gamma)./(1+(b01.*P.*exp(delU1./(8.314.*T)))^gamma) ...
                            + qs2.*((b02.*P.*exp(delU2./(8.314.*T)))^gamma)./(1+(b02.*P.*exp(delU2./(8.314.*T)))^gamma);
                    case 'TOTH'
                        qvals(jj,kk) = qs1.*b01.*P.*exp(delU1./(8.314.*T))/(1+(b01.*P.*exp(delU1./(8.314.*T))).^toth).^(1./toth);
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
            xlim([0 max(x)])
            ylim([0 1.1.*max(z)])
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
    case {'VIRIAL','VIRIAL2'}
        isothermData.isothermFit = [headerRow;qvals(1,:)' lnPvals];
        isothermData.experiment = [z log(x) y];
        isothermData.heatAdsorption = [qvals' delH'];
    otherwise
        if flagConcUnits
            isothermData.isothermFit = [];
            isothermData.isothermFit = [linspace(0,x(find(x==max(max(x))))./(1e5./(8.314.*y(find(x==max(max(x)))))),length(qvals))' qvals];
            isothermData.confidenceBounds = [uncBounds(:,1)./(1e5./(8.314.*uncBounds(find(x==max(max(x))),4))) uncBounds(:,2) uncBounds(:,3) uncBounds(:,4)];
            outScatter(:,1) = outScatter(:,1)./(1e5./(8.314.*outScatter(:,3)));
            isothermData.confidenceRegion = outScatter;
        else
            isothermData.isothermFit = [headerRow;Pvals(1,:)' qvals];
            isothermData.confidenceRegion = outScatter;
            isothermData.confidenceBounds = uncBounds;
        end
end
isothermData.isothermParameters = [parsDisp' conRange95Disp];
isothermData.gitCommitID = gitCommitID;

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
    switch isothermModel
        case {'VIRIAL','VIRIAL2'}
        otherwise
            figure
            plot(isothermData.isothermFit(:,1),isothermData.isothermFit(:,2:end),'-k','LineWidth',1.5)
            hold on;
            plot(isothermData.experiment(:,1),isothermData.experiment(:,2),'ok')
            scatter(isothermData.confidenceBounds(:,1),isothermData.confidenceBounds(:,2),0.5,'MarkerEdgeColor','b','MarkerEdgeAlpha',0.5);
            scatter(isothermData.confidenceBounds(:,1),isothermData.confidenceBounds(:,3),0.5,'MarkerEdgeColor','b','MarkerEdgeAlpha',0.5);
            xlabel('Pressure [bar]');
            ylabel('Adsorbed amount [mol/kg]');
            xlim([0 x(find(x==max(max(x))))./(1e5./(8.314.*y(find(x==max(max(x))))))])
            ylim([0 1.1.*max(z)])
            box on
            set(gca,'YScale','linear','XScale','linear','FontSize',15,'LineWidth',1)
            grid on; axis square
    end
end

varargout{1} = isothermData;
if strcmp(currentDir(end),'ERASE')
    cd ..
end
end