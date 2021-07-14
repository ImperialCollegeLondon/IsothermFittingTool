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
% Generating objective function for maximum log-likelihood estimator for
% fitting isotherm parameters
%
% Last modified:
% - 2021-03-11, HA: Initial creation
%
% Input arguments:
% - x, y, z:             Pressure (x), Temperature (y), and adsorbed amount
%                        (z) data from experiments
%
% - nbins:               Number of groups data is binned to based on
%                        pressure (x) range (1 for MLE)
%
% - isothermModel:       Isotherm model for which the data needs to be
%                        fitted ('DSL' or 'DSS')
%
% - varargin:            Arguments including the isotherm parameter
%                        variables for prediction
%
% Output arguments:
% - objfun:              Objective function for optimisation via MLE method
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function objfun = generateMLEfun(x,y,z, nbins, isothermModel, isoRef, varargin)
% Discretize the data range into 'nbins' evenly spaced bins of pressure
% ranges ('nbins = 1' for MLE')
[bins] = discretize(x,nbins);
% Create vector for storing sum of error for each bin
err = zeros(nbins,1);
% Number of data points
Nt = length(bins);

expData = [x,z,y];
expData = sortrows(expData,3);

x = expData(:,1);
z = expData(:,2);
y = expData(:,3);

temperatureValues = unique(y);
qRefIndexTemp = zeros(length(temperatureValues),1);
for ii = 1:length(temperatureValues)
    qRefIndexTemp(ii,1) = find(y == temperatureValues(ii),1,'first');
    qRefIndexTemp(ii,2) = find(y == temperatureValues(ii),1,'last');
end

% Find qref for the experimental data
qRefMax = max(z(qRefIndexTemp(:,2)));
qRefTemp = z(qRefIndexTemp(:,2));
normalizationFactorTemp = qRefMax./qRefTemp;
normalizationFactor = zeros(length(x),1);

for ii = 1:length(temperatureValues)
    normalizationFactor(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2),1) = normalizationFactorTemp(ii);
end

% Calculate error based on isotherm model
switch isothermModel
    % Calculate error for DSL model
    case 'DSL'
        qs1 = varargin{1}.*isoRef(1);
        qs2 = varargin{2}.*isoRef(2);
        b01 = varargin{3}.*isoRef(3);
        b02 = varargin{4}.*isoRef(4);
        delU1 = varargin{5}.*isoRef(5);
        delU2 = varargin{6}.*isoRef(6);
        
        % Loop for calculating the sum of errors for each bin
        for jj = 1:length(err)
            for kk = 1:length(bins)
                qfun = qs1*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                    + qs2*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))));
                if bins(kk) == jj
                    err(jj) = (err(jj) + (normalizationFactor(kk).*(z(kk) - qfun))^2);
                end
            end
            err(jj) = err(jj);
        end
    % Calculate error for DSS model
    case 'DSS'
        qs1 = varargin{1}.*isoRef(1);
        qs2 = varargin{2}.*isoRef(2);
        b01 = varargin{3}.*isoRef(3);
        b02 = varargin{4}.*isoRef(4);
        delU1 = varargin{5}.*isoRef(5);
        delU2 = varargin{6}.*isoRef(6);
        gamma = varargin{7}.*isoRef(7);

        % Loop for calculating the sum of errors for each bin
        for jj = 1:length(err)
            for kk = 1:length(bins)
                qfun = qs1*(b01*x(kk)*exp(delU1/(8.314*y(kk))))^gamma/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^gamma) ...
                    + qs2*(b02*x(kk)*exp(delU2/(8.314*y(kk))))^gamma/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk))))^gamma);
                if bins(kk) == jj
                    err(jj) = (err(jj) + (normalizationFactor(kk).*(z(kk) - qfun))^2);
                end
            end
            err(jj) = err(jj);
        end
end
% Define the objective function as the output
objfun = (Nt/2 * log(sum(err)));
end