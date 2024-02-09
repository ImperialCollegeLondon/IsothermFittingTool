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
% Calculate error based on isotherm model
switch isothermModel
    % Calculate error for Universal isotherm model for Type VI
    case 'UNIV6'
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
        % [qs a1 a2 a3 e01 e02 e03 e04 m1 m2 m3 m4]
        qs =  varargin{1}.*isoRef(1);
        a1 =  varargin{2}.*isoRef(2);
        a2 =  varargin{3}.*isoRef(3);
        a3 =  varargin{4}.*isoRef(4);
        e01 = varargin{5}.*isoRef(5);
        e02 = varargin{6}.*isoRef(6);
        e03 = varargin{7}.*isoRef(7);
        e04 = varargin{8}.*isoRef(8);
        m1 =  varargin{9}.*isoRef(9);
        m2 = varargin{10}.*isoRef(10);
        m3 = varargin{11}.*isoRef(11);
        m4 = varargin{12}.*isoRef(12);
        ps = varargin{13}.*isoRef(13);

        for jj = 1:length(err)
            for kk = 1:length(bins)
                if a1 + a2 + a3 > 1
                    qfun = 0;
                else
                    qfun = computeUNIV6Loading(x(kk),y(kk), qs, a1, a2, a3, e01, e02, e03, e04, m1, m2, m3, m4, ps);
                end
                if bins(kk) == jj
                    err(jj) = (err(jj) + (normalizationFactor(kk).*(z(kk) - qfun))^2);
                end
            end
        end
        err;
        % Calculate error for Statistical isotherm model for zeolites
    case 'STATZ'
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

        omega = round(varargin{1}).*isoRef(1);
        beta = varargin{2}.*isoRef(2);
        b01 = varargin{3}.*isoRef(3);
        delU1 = varargin{4}.*isoRef(4);
        vc = varargin{5};
        vm = varargin{6};
        Na = 6.022e20; % Avogadros constant [molecules/mmol]
        for jj = 1:length(err)
            for kk = 1:length(bins)
                if omega > vc/beta
                    qfun = 0;
                else
                    qfun = computeStatZLoading(x(kk),y(kk),b01,delU1,beta,omega,vc);
                end
                if bins(kk) == jj
                    err(jj) = (err(jj) + (normalizationFactor(kk).*(vm./(vc.*Na)).*(z(kk) - qfun))^2);
                end
            end
        end
        % Calculate error for Statistical isotherm model for zeolites
    case 'STATZGATE'
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
        Na = 6.022e20; % Avogadros constant [molecules/mmol]
        omega = varargin{1}.*isoRef(1);
        beta = varargin{2}.*isoRef(2);
        b01 = varargin{3}.*isoRef(3);
        delU1 = varargin{4}.*isoRef(4);
        kgate = varargin{5}.*isoRef(5);
        cgate = varargin{6}.*isoRef(6);
        gamma = varargin{7}.*isoRef(7);
        vc1 = varargin{8};
        vc2 = varargin{9};
        vm  = varargin{10};
        for jj = 1:length(err)
            for kk = 1:length(bins)
                if omega > vc2/beta
                    qfun = 0;
                else
                    qfun = computeStatZGATELoading(x(kk),y(kk),b01,delU1,beta,omega,kgate,cgate,gamma,vc1,vc2);
                end
                if bins(kk) == jj
                    err(jj) = (err(jj) + (normalizationFactor(kk).*(vm./(vc2.*Na)).*(z(kk) - qfun))^2);
                end
            end
        end

        % Calculate error for DSL model
    case 'DSL'
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
    case 'DSL2'
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

        qs1a = varargin{1}.*isoRef(1);
        qs2a = varargin{2}.*isoRef(2);
        qs1b = varargin{3}.*isoRef(3);
        qs2b = varargin{4}.*isoRef(4);
        b01 = varargin{5}.*isoRef(5);
        b02 = varargin{6}.*isoRef(6);
        delU1 = varargin{7}.*isoRef(7);
        delU2 = varargin{8}.*isoRef(8);

        % Loop for calculating the sum of errors for each bin
        for jj = 1:length(err)
            for kk = 1:length(bins)
                qs1 = qs1a + qs1b./y(kk);
                qs2 = qs2a + qs2b./y(kk);
                qfun = qs1*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                    + qs2*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))));
                if bins(kk) == jj
                    err(jj) = (err(jj) + (normalizationFactor(kk).*(z(kk) - qfun))^2);
                end
            end
            err(jj) = err(jj);
        end
        % Calculate error for DSL model
    case 'HDSL'
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

        qs1   = varargin{1}.*isoRef(1);
        qs2   = varargin{2}.*isoRef(2);
        b01   = varargin{3}.*isoRef(3);
        b02   = varargin{4}.*isoRef(4);
        delU1 = varargin{5}.*isoRef(5);
        delU2 = varargin{6}.*isoRef(6);
        qsH   = varargin{7}.*isoRef(7);
        b0H   = varargin{8}.*isoRef(8);
        delUH = varargin{9}.*isoRef(9);

        % Loop for calculating the sum of errors for each bin
        for jj = 1:length(err)
            for kk = 1:length(bins)
                qfun = qs1*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                    + qs2*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk))))) ...
                    + qsH*(b0H*x(kk)*exp(delUH/(8.314*y(kk))));
                if bins(kk) == jj
                    err(jj) = (err(jj) + (normalizationFactor(kk).*(z(kk) - qfun))^2);
                end
            end
            err(jj) = err(jj);
        end
    case 'TSL'
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

        qs1   = varargin{1}.*isoRef(1);
        qs2   = varargin{2}.*isoRef(2);
        b01   = varargin{3}.*isoRef(3);
        b02   = varargin{4}.*isoRef(4);
        delU1 = varargin{5}.*isoRef(5);
        delU2 = varargin{6}.*isoRef(6);
        qs3   = varargin{7}.*isoRef(7);
        b03   = varargin{8}.*isoRef(8);
        delU3 = varargin{9}.*isoRef(9);

        % Loop for calculating the sum of errors for each bin
        for jj = 1:length(err)
            for kk = 1:length(bins)
                qfun = qs1*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                    + qs2*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk))))) ...
                    + qs3*(b03*x(kk)*exp(delU3/(8.314*y(kk))))/(1+(b03*x(kk)*exp(delU3/(8.314*y(kk)))));
                if bins(kk) == jj
                    err(jj) = (err(jj) + (normalizationFactor(kk).*(z(kk) - qfun))^2);
                end
            end
            err(jj) = err(jj);
        end
        % Calculate error for DSS model
    case 'DSS'
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
    case 'DSS2'
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
                qfun = qs1*(b01*x(kk)^gamma*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)^gamma*exp(delU1/(8.314*y(kk))))) ...
                    + qs2*(b02*x(kk)^gamma*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)^gamma*exp(delU2/(8.314*y(kk)))));
                if bins(kk) == jj
                    err(jj) = (err(jj) + (normalizationFactor(kk).*(z(kk) - qfun))^2);
                end
            end
            err(jj) = err(jj);
        end
    case 'TOTH'
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

        qs1 = varargin{1}.*isoRef(1);
        b01 = varargin{3}.*isoRef(3);
        delU1 = varargin{5}.*isoRef(5);
        toth = varargin{7}.*isoRef(7);

        % Loop for calculating the sum of errors for each bin
        for jj = 1:length(err)
            for kk = 1:length(bins)
                qfun = qs1.*b01.*x(kk).*exp(delU1./(8.314.*y(kk)))/(1+(b01.*x(kk).*exp(delU1./(8.314.*y(kk)))).^toth).^(1./toth);
                if bins(kk) == jj
                    err(jj) = (err(jj) + (normalizationFactor(kk).*(z(kk) - qfun))^2);
                end
            end
            err(jj) = err(jj);
        end
    case 'TOTH2'
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

        qs1 = varargin{1}.*isoRef(1);
        b01 = varargin{3}.*isoRef(3);
        delU1 = varargin{5}.*isoRef(5);
        toth0 = varargin{7}.*isoRef(7);
        totha = varargin{8}.*isoRef(8);

        % Loop for calculating the sum of errors for each bin
        for jj = 1:length(err)
            for kk = 1:length(bins)
                toth = toth0 + totha.*(1-298.15./y(kk));
                qfun = qs1.*b01.*x(kk).*exp(delU1./(8.314.*y(kk)))/(1+(b01.*x(kk).*exp(delU1./(8.314.*y(kk)))).^toth).^(1./toth);
                if bins(kk) == jj
                    err(jj) = (err(jj) + (normalizationFactor(kk).*(z(kk) - qfun))^2);
                end
            end
            err(jj) = err(jj);
        end
    case 'TOTH3'
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

        qs10 = varargin{1}.*isoRef(1);
        b01 = varargin{3}.*isoRef(3);
        delU1 = varargin{5}.*isoRef(5);
        toth0 = varargin{7}.*isoRef(7);
        totha = varargin{8}.*isoRef(8);
        chi = varargin{9}.*isoRef(9);

        % Loop for calculating the sum of errors for each bin
        for jj = 1:length(err)
            for kk = 1:length(bins)
                toth = toth0 + totha.*(1-298.15./y(kk));
                qs1 = qs10.*exp(chi*(1-y(kk)./298.15));
                qfun = qs1.*b01.*x(kk).*exp(delU1./(8.314.*y(kk)))/(1+(b01.*x(kk).*exp(delU1./(8.314.*y(kk)))).^toth).^(1./toth);
                if bins(kk) == jj
                    err(jj) = (err(jj) + (normalizationFactor(kk).*(z(kk) - qfun))^2);
                end
            end
            err(jj) = err(jj);
        end
    case 'TOTHCHEM'
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

        qs10 = varargin{1}.*isoRef(1);
        b01 = exp(varargin{3}.*isoRef(3));
        delU1 = varargin{5}.*isoRef(5);
        toth0 = varargin{7}.*isoRef(7);
        totha = varargin{8}.*isoRef(8);
        chi = varargin{9}.*isoRef(9);
        qsC = varargin{10}.*isoRef(10);
        b0C = exp(varargin{11}.*isoRef(11));
        delUC = varargin{12}.*isoRef(12);
        EaC = varargin{13}.*isoRef(13);

        % Loop for calculating the sum of errors for each bin
        for jj = 1:length(err)
            for kk = 1:length(bins)
                toth = toth0 + totha.*(1-298.15./y(kk));
                qs1 = qs10.*exp(chi*(1-y(kk)./298.15));
                qfun = qs1.*b01.*x(kk).*exp(delU1./(8.314.*y(kk)))/(1+(b01.*x(kk).*exp(delU1./(8.314.*y(kk)))).^toth).^(1./toth) ...
                    + exp(-EaC./(8.314.*y(kk))).*qsC.*b0C.*x(kk).*exp(delUC./(8.314.*y(kk)))./(1+b0C.*x(kk).*exp(delUC./(8.314.*y(kk))));
                if bins(kk) == jj
                    err(jj) = (err(jj) + (normalizationFactor(kk).*(z(kk) - qfun))^2);
                end
            end
            err(jj) = err(jj);
        end
    case 'GAB'
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

        qs1  = varargin{1}.*isoRef(1);
        parC = varargin{2}.*isoRef(2);
        parD = varargin{3}.*isoRef(3);
        parF = varargin{4}.*isoRef(4);
        parG = varargin{5}.*isoRef(5);

        % Loop for calculating the sum of errors for each bin
        for jj = 1:length(err)
            for kk = 1:length(bins)
                qfun = computeGABLoading(x(kk),y(kk),qs1,parC,parD,parF,parG);
                if bins(kk) == jj
                    err(jj) = (err(jj) + (normalizationFactor(kk).*(z(kk) - qfun))^2);
                end
            end
            err(jj) = err(jj);
        end
        % for Virial model
    case 'VIRIAL'
        expData = [x,z,y];
        expData = sortrows(expData,1);
        expData = expData(length(unique(y))+1:end,:);
        expData = sortrows(expData,3);

        x = expData(:,1);
        z = expData(:,2);
        y = expData(:,3);

        temperatureValues = unique(y);
        lnPRefIndexTemp = zeros(length(temperatureValues),1);
        for ii = 1:length(temperatureValues)
            lnPRefIndexTemp(ii,1) = find(y == temperatureValues(ii),1,'first');
            lnPRefIndexTemp(ii,2) = find(y == temperatureValues(ii),1,'last');
        end

        % Find qref for the experimental data
        lnPRefMax = max(log(x(lnPRefIndexTemp(:,2))));
        lnPRefTemp = log(x(lnPRefIndexTemp(:,2)));
        normalizationFactorTemp = lnPRefMax./lnPRefTemp;
        normalizationFactor = zeros(length(x),1);

        for ii = 1:length(temperatureValues)
            normalizationFactor(lnPRefIndexTemp(ii,1):lnPRefIndexTemp(ii,2),1) = normalizationFactorTemp(ii);
        end

        % Discretize the data range into 'nbins' evenly spaced bins of pressure
        % ranges ('nbins = 1' for MLE')
        [bins] = discretize(x,nbins);
        % Create vector for storing sum of error for each bin
        err = zeros(nbins,1);
        % Number of data points
        Nt = length(bins);

        a0 = varargin{1}.*isoRef(1);
        a1 = varargin{2}.*isoRef(2);
        a2 = varargin{3}.*isoRef(3);
        a3 = varargin{4}.*isoRef(4);
        b0 = varargin{5}.*isoRef(5);
        b1 = varargin{6}.*isoRef(6);
        b2 = 0;
        b3 = 0;

        parameters = [a0, a1, a2, a3, b0, b1, b2, b3];

        % Loop for calculating the sum of errors for each bin
        for jj = 1:length(err)
            for kk = 1:length(bins)
                lnP = log(z(kk)) + 1/y(kk).*(parameters(1) + parameters(2)*z(kk) + ...
                    parameters(3)*z(kk)^2 + parameters(4)*z(kk)^3) + parameters(5) ...
                    + parameters(6)*z(kk)+ parameters(7)*z(kk).^2+ parameters(8)*z(kk).^3;
                if bins(kk) == jj
                    err(jj) = err(jj) + (normalizationFactor(kk).*(log(x(kk)) - lnP))^2;
                end
            end
            err(jj) = err(jj);
        end
    case 'VIRIAL2'
        expData = [x,z,y];
        expData = sortrows(expData,1);
        expData = expData(length(unique(y))+1:end,:);
        expData = sortrows(expData,3);

        x = expData(:,1);
        z = expData(:,2);
        y = expData(:,3);

        temperatureValues = unique(y);
        lnPRefIndexTemp = zeros(length(temperatureValues),1);
        for ii = 1:length(temperatureValues)
            lnPRefIndexTemp(ii,1) = find(y == temperatureValues(ii),1,'first');
            lnPRefIndexTemp(ii,2) = find(y == temperatureValues(ii),1,'last');
        end

        % Find qref for the experimental data
        lnPRefMax = max(log(x(lnPRefIndexTemp(:,2))));
        lnPRefTemp = log(x(lnPRefIndexTemp(:,2)));
        normalizationFactorTemp = lnPRefMax./lnPRefTemp;
        normalizationFactor = zeros(length(x),1);

        for ii = 1:length(temperatureValues)
            normalizationFactor(lnPRefIndexTemp(ii,1):lnPRefIndexTemp(ii,2),1) = normalizationFactorTemp(ii);
        end

        % Discretize the data range into 'nbins' evenly spaced bins of pressure
        % ranges ('nbins = 1' for MLE')
        [bins] = discretize(x,nbins);
        % Create vector for storing sum of error for each bin
        err = zeros(nbins,1);
        % Number of data points
        Nt = length(bins);

        a0 = varargin{1}.*isoRef(1);
        a1 = varargin{2}.*isoRef(2);
        a2 = varargin{3}.*isoRef(3);
        a3 = varargin{4}.*isoRef(4);
        b0 = varargin{5}.*isoRef(5);
        b1 = varargin{6}.*isoRef(6);
        b2 = varargin{7}.*isoRef(7);
        b3 = varargin{8}.*isoRef(8);

        parameters = [a0, a1, a2, a3, b0, b1, b2, b3];

        % Loop for calculating the sum of errors for each bin
        for jj = 1:length(err)
            for kk = 1:length(bins)
                lnP = log(z(kk)) + 1/y(kk).*(parameters(1) + parameters(2)*z(kk) + ...
                    parameters(3)*z(kk)^2 + parameters(4)*z(kk)^3) + parameters(5) ...
                    + parameters(6)*z(kk)+ parameters(7)*z(kk).^2+ parameters(8)*z(kk).^3;
                if bins(kk) == jj
                    err(jj) = err(jj) + (normalizationFactor(kk).*(log(x(kk)) - lnP))^2;
                end
            end
            err(jj) = err(jj);
        end
end
% Define the objective function as the output
objfun = (Nt/2 * log(sum(err)));
end