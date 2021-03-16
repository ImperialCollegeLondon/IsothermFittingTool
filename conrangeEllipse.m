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
% Determining the 95% confidence regions for parameters estimated using MLE
% or WSS for a data set assumed to follow a normal distribution
%
% Last modified:
% - 2021-03-15, HA: Updated method to obtain confidence intervals and made
% corrections to the code
% - 2021-03-11, HA: Initial creation
%
% Input arguments:
% - x, y, z:             Pressure (x), Temperature (y), and adsorbed amount
%                        (z) data from experiments
%
% - qfit:                Adsorbed amount calculated using fitted parameters
%                        for either DSS or DSL model
%
% - isothermModel:       Isotherm model for which the data needs to be
%                        fitted ('DSL' or 'DSS')
%
% - varargin:            Arguments including the isotherm parameter
%                        variables for prediction
%
% Output arguments:
% - conRange95:          Vector of confidence ranges for each parameter
%
% - fval:                Value of the objective function at the optimal
%                        point
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [conRange95] = conrangeEllipse(x,y,z,  qfit, isothermModel, varargin)
% Calculate standard deviation of the data (not needed)
stDevData = 1/length(x) * sum((z-qfit).^2);
% Generate and solve global optimisation problem for confidence regions
% based on isotherm model
switch isothermModel
    % For DSL model
    case 'DSL'
        qs1 = varargin{1};
        qs2 = varargin{2};
        b01 = varargin{3};
        b02 = varargin{4};
        delU1 = varargin{5};
        delU2 = varargin{6};
        % Create empty sensitivity matrix
        sensitivityMatrix = zeros(length(x),6);
        % degree of variation of parameter to calculate sensitivity
        del = 0.0001;
        % Calculate sensitivity at every data point for each parameter
        for jj = 1:6
            for kk = 1:length(x)
                if jj == 1
                    sensitivityMatrix(kk,jj) = ((((1+del)*qs1)*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                        + qs2*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))))) - qfit(kk))/(del*qs1);
                elseif jj == 2
                    sensitivityMatrix(kk,jj) = ((qs1*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                        + ((1+del)*qs2)*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))))) - qfit(kk))/(del*qs2);
                elseif jj == 3
                    sensitivityMatrix(kk,jj) = ((qs1*(((1+del)*b01)*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(((1+del)*b01)*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                        + qs2*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))))) - qfit(kk))/(del*b01);
                elseif jj == 4
                    sensitivityMatrix(kk,jj) = ((qs1*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                        + qs2*(((1+del)*b02)*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(((1+del)*b02)*x(kk)*exp(delU2/(8.314*y(kk)))))) - qfit(kk))/(del*b02);
                elseif jj == 5
                    sensitivityMatrix(kk,jj) = ((qs1*(b01*x(kk)*exp((1+del)*delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp((1+del)*delU1/(8.314*y(kk))))) ...
                        + qs2*((b02)*x(kk)*exp(delU2/(8.314*y(kk))))/(1+((b02)*x(kk)*exp(delU2/(8.314*y(kk)))))) - qfit(kk))/(del*delU1);
                else
                    sensitivityMatrix(kk,jj) = ((qs1*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                        + qs2*((b02)*x(kk)*exp((1+del)*delU2/(8.314*y(kk))))/(1+((b02)*x(kk)*exp((1+del)*delU2/(8.314*y(kk)))))) - qfit(kk))/(del*delU2);
                end
            end
        end
        % estimated Hessian Matrix for the data set (Non-linear parameter estimation
        % by Yonathan Bard (1974) pg. 178
        hessianMatrix = 1/stDevData*transpose(sensitivityMatrix)*sensitivityMatrix;
        % Confidence range given by chi squared distribution at Np degrees
        % of freedom (independent parameter conf intervals)
        conRange95 = sqrt(chi2inv(0.95,6)./diag(hessianMatrix));
        
        % for DSS model
    case 'DSS'
        qs1 = varargin{1};
        qs2 = varargin{2};
        b01 = varargin{3};
        b02 = varargin{4};
        delU1 = varargin{5};
        delU2 = varargin{6};
        gamma = varargin{7};
        % Create empty sensitivity matrix
        sensitivityMatrix = zeros(length(x),7);
        % degree of variation of parameter to calculate sensitivity
        del = 0.0001;
        % Calculate sensitivity at every data point for each parameter
        for jj = 1:7
            for kk = 1:length(x)
                if jj == 1
                    sensitivityMatrix(kk,jj) = ((((1+del)*qs1)*(b01*x(kk)*exp(delU1/(8.314*y(kk))))^gamma/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^gamma) ...
                        + qs2*(b02*x(kk)*exp(delU2/(8.314*y(kk))))^gamma/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk))))^gamma)) - qfit(kk))/(del*qs1);
                elseif jj == 2
                    sensitivityMatrix(kk,jj) = ((qs1*(b01*x(kk)*exp(delU1/(8.314*y(kk))))^gamma/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^gamma) ...
                        + ((1+del)*qs2)*(b02*x(kk)*exp(delU2/(8.314*y(kk))))^gamma/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk))))^gamma)) - qfit(kk))/(del*qs2);
                elseif jj == 3
                    sensitivityMatrix(kk,jj) = ((qs1*(((1+del)*b01)*x(kk)*exp(delU1/(8.314*y(kk))))^gamma/(1+(((1+del)*b01)*x(kk)*exp(delU1/(8.314*y(kk))))^gamma) ...
                        + qs2*(b02*x(kk)*exp(delU2/(8.314*y(kk))))^gamma/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk))))^gamma)) - qfit(kk))/(del*b01);
                elseif jj == 4
                    sensitivityMatrix(kk,jj) = ((qs1*(b01*x(kk)*exp(delU1/(8.314*y(kk))))^gamma/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^gamma) ...
                        + qs2*(((1+del)*b02)*x(kk)*exp(delU2/(8.314*y(kk))))^gamma/(1+(((1+del)*b02)*x(kk)*exp(delU2/(8.314*y(kk))))^gamma)) - qfit(kk))/(del*b02);
                elseif jj == 5
                    sensitivityMatrix(kk,jj) = ((qs1*(b01*x(kk)*exp((1+del)*delU1/(8.314*y(kk))))^gamma/(1+(b01*x(kk)*exp((1+del)*delU1/(8.314*y(kk))))^gamma) ...
                        + qs2*((b02)*x(kk)*exp(delU2/(8.314*y(kk))))^gamma/(1+((b02)*x(kk)*exp(delU2/(8.314*y(kk))))^gamma)) - qfit(kk))/(del*delU1);
                elseif jj == 6
                    sensitivityMatrix(kk,jj) = ((qs1*(b01*x(kk)*exp(delU1/(8.314*y(kk))))^gamma/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^gamma) ...
                        + qs2*((b02)*x(kk)*exp((1+del)*delU2/(8.314*y(kk))))^gamma/(1+((b02)*x(kk)*exp((1+del)*delU2/(8.314*y(kk))))^gamma)) - qfit(kk))/(del*delU2);
                else
                    sensitivityMatrix(kk,jj) = ((qs1*(b01*x(kk)*exp(delU1/(8.314*y(kk))))^((1+del)*gamma)/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^((1+del)*gamma)) ...
                        + qs2*((b02)*x(kk)*exp(delU2/(8.314*y(kk))))^((1+del)*gamma)/(1+((b02)*x(kk)*exp(delU2/(8.314*y(kk))))^((1+del)*gamma))) - qfit(kk))/(del*gamma);
                end
            end
        end
        % Hessian Matrix for the data set (Non-linear parameter estimation
        % by Yonathan Bard (1974) pg. 178)
        hessianMatrix = 1/stDevData*transpose(sensitivityMatrix)*sensitivityMatrix;
        % Confidence range given by chi squared distribution at Np degrees
        % of freedom (independent parameter conf intervals)
        conRange95 = sqrt(chi2inv(0.95,7)./diag(hessianMatrix));
end
end