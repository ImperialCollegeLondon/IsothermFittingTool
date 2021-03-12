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
% Fits input experimental adsorption isotherm data to either dual-site
% Langmuir or dual-site Sips adsorption model, and outputs isotherm
% parameters and ellipsoidal confidence bounds
%
% Last modified:
% - 2021-03-11, HA: Initial creation
%
% Input arguments:
%
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALISATION
% Clear command window and workspace
% clc
clear
% Load input experimental data from *.mat or *.csv file in a 3 column
% format with Pressure (bar), adsorbed amount (-), temperature (K) in the 3
% columns respectively
load zif8Data
fitData = zif8Data;
% Determine number of bins you want the experimental data to be binned to
% in terms of the total pressure range
nbins = 3;
% Pressure (x), adsorbed amount (z), Temperature (y) data from input
x = fitData(:,1);
z = fitData(:,2);
y = fitData(:,3);
% Select isotherm model for fitting
isothermModel = 'DSS'; % DSL = Dual site Langmuir. DSS = Dual site Sips
% Select fitting method
fittingMethod = 'WSS'; % WSS = weighted sum of squares, MLE = max log likelihood estimator

%% GENERATING AND SOLVING OPTIMISATION PROBLEM
% Set fitting procedure based on isotherm model
switch isothermModel
    case 'DSL'
        rng default % For reproducibility
        % Create gs, a GlobalSearch solver with its properties set to the defaults.
        gs = GlobalSearch('NumTrialPoints',1400,'NumStageOnePoints',200);
        % Set objective function based on fitting method
        switch fittingMethod
            case 'MLE'
                % Number of bins is automatically set to 1 for MLE as
                % weights cannot be assigned in MLE
                nbins =1;
                % Generate objective function for MLE method
                optfunc = @(par) generateMLEfun(x, y, z, nbins, isothermModel, par(1), par(2), par(3), ...
                          par(4), par(5), par(6));
            case 'WSS'
                % Generate objective function for WSS method
                optfunc = @(par) generateWSSfun(x, y, z, nbins, isothermModel, par(1), par(2), par(3), ...
                          par(4), par(5), par(6));
        end
        % Initial conditions, lower bounds, and upper bounds for parameters
        % in DSL isotherm model
        x0 = [3,3,0.001,0.001,2e4,2e4];
        lb = [0,0,0,0,0,0];
        ub = [20,20,1,1,8e4,8e4];
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
        % Calculate fitted isotherm loadings for conditions (P,T)
        % corresponding to experimental data
        qfit  = qs1.*(b01.*x.*exp(delU1./(8.314.*y)))./(1+(b01.*x.*exp(delU1./(8.314.*y)))) ...
            + qs2.*(b02.*x.*exp(delU2./(8.314.*y)))./(1+(b02.*x.*exp(delU2./(8.314.*y))));
        % Calculate ellipsoidal confidence regions (delta parameter) for
        % fitted parameters
        [conRange95, fvalConf] = conrangeEllipse(x, y, z, qfit, isothermModel, qs1, qs2, b01, b02, delU1, delU2);
        % Convert confidence regions to percentage error
        percentageError = conRange95./parVals *100
        
    case 'DSS'
        rng default % For reproducibility
        % Create gs, a GlobalSearch solver with its properties set to the defaults.
        gs = GlobalSearch('NumTrialPoints',1400,'NumStageOnePoints',200);
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
        x0 = [3,3,0.001,0.001,2e4,2e4,1];
        lb = [0,0,0,0,0,0,0];
        ub = [20,20,1,1,+inf,+inf,2];
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
        % Calculate ellipsoidal confidence regions (delta parameter) for
        % fitted parameters
        [conRange95, fvalConf] = conrangeEllipse(x, y, z, qfit, isothermModel, qs1, qs2, b01, b02, delU1, delU2, gamma);
        % Convert confidence regions to percentage error
        percentageError = conRange95./parVals *100
end

%% PLOT RESULTING DATA
figure
% plot of experimental data and fitted data (q vs P)
subplot(1,3,1)
semilogx(x,z,'ok');
hold on
semilogx(x,qfit,'ob');
xlabel('Pressure [bar]');
ylabel('q [mol/kg]');
% quantile-quantile plot of experimental data vs normal distribution
subplot(1,3,2)
pd = fitdist(z-qfit,'Normal');
qqplot(qfit,pd);
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

Z = generateObjfunContour(x,y,z,nbins,isothermModel,fittingMethod,parVals);