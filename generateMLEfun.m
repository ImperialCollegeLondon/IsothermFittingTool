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
        Na = 6.022e20; % Avogadros constant [molecules/mmol]
        % Number of data points
        Nt = length(bins);
        omega = round(varargin{1}).*isoRef(1);
        beta = varargin{2}.*isoRef(2);
        b01 = varargin{3}.*isoRef(3);
        delU1 = varargin{4}.*isoRef(4);
        vc = varargin{5};
        vm = varargin{6};
        % cutoffq = varargin{7};
        % cutoffq1 = cutoffq(1).*((vc.*Na)./(vm));
        % cutoffq2 = cutoffq(2).*((vc.*Na)./(vm));
        % cutoffq0 = cutoffq(3).*((vc.*Na)./(vm));

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
        % normalizationFactorTemp = zeros(length(zeros(length(x),1)),1);
        % normalizationFactorTemp(:) = 1;
        normalizationFactor = zeros(length(x),1);

        for ii = 1:length(temperatureValues)
            % normalizationFactor(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2),1) = normalizationFactorTemp(ii)./length(normalizationFactor(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2),1)).*length(x);
            normalizationFactor(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2),1) = normalizationFactorTemp(ii);
            % normalizationFactor(qRefIndexTemp(ii,1):qRefIndexTemp(ii,1)+5,1) = normalizationFactorTemp(ii);
        end
        
        % indexCutoff = zeros(length(temperatureValues),2);

        % if cutoffq0 > 0
            % for ii = 1:length(temperatureValues)
            %     % indexCutoff(ii,1) = qRefIndexTemp(ii,1);
            %     % if ~isempty(find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) < cutoffq0,1,"last")) & ~isempty(find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) > cutoffq1,1,"first"))
            %     %     indexCutoff(ii,2) = indexCutoff(ii,1) + find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) < cutoffq0,1,"last");
            %     %     normalizationFactor(indexCutoff(ii,1):indexCutoff(ii,2)-1) = 1.*z(indexCutoff(ii,2)+1)./z(indexCutoff(ii,2)-1).*normalizationFactor(indexCutoff(ii,1):indexCutoff(ii,2)-1);
            %     % else
            %     %     % indexCutoff(ii,2) = indexCutoff(ii,1) + find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) < cutoffq0,1,"last");
            %     %     normalizationFactor(indexCutoff(ii,1):indexCutoff(ii,2)-1) = 1.*z(indexCutoff(ii-1,2)+1)./z(indexCutoff(ii-1,2)-1).*normalizationFactor(indexCutoff(ii,1):indexCutoff(ii,2)-1);
            %     % end
            %     % normalizationFactor(indexCutoff(ii,1):indexCutoff(ii,2),1)     = max(z_new(indexCutoff(ii,1):indexCutoff(ii,2),1))./z_new(indexCutoff(ii,1):indexCutoff(ii,2),1);
            %     % normalizationFactor(indexCutoff(ii,2)+1:qRefIndexTemp(ii,2),1) = ((z_new(indexCutoff(ii,2)+1:qRefIndexTemp(ii,2))))./max(z_new(indexCutoff(ii,2)+1:qRefIndexTemp(ii,2),1));
            %     % normalizationFactor(indexCutoff(ii,1):indexCutoff(ii,2),1) = z_new(indexCutoff(ii,1):indexCutoff(ii,2))./cutoffq0;
                % % normalizationFactor(indexCutoff(ii,1):indexCutoff(ii,2)-1) = z(indexCutoff(ii,2)+1)./z(indexCutoff(ii,2)-1).*normalizationFactor(indexCutoff(ii,1):indexCutoff(ii,2)-1);
            % end
        % else
            % normalizationFactor(:) = 1;
        % end

        % figure
        % scatter(x,z.*normalizationFactor)
        % hold on
        % scatter(x,z)

        
        for jj = 1:length(err)
            for kk = 1:length(bins)
                % if omega > vc/beta
                %     qfun = 0;
                % else
                    qfun = computeStatZLoading(x(kk),y(kk),b01,delU1,beta,omega,vc);
                % end
                if bins(kk) == jj
                    % if kk < round(length(err)./10)
                        % err(jj) = (err(jj) + (normalizationFactor(kk).*(vm./(vc.*Na)).*(z(kk) - qfun))^2);
                        % err(jj) = (err(jj) + (normalizationFactor(kk).*(z(kk)./((vc.*Na)./(vm)) - qfun./((vc.*Na)./(vm))))^2);
                        err(jj) = (err(jj) + (normalizationFactor(kk).*(z(kk) - qfun))^2);
                    % else
                        % err(jj) = (err(jj) + (normalizationFactor(kk).*(vm./(vc.*Na)).*(z(kk) - qfun))^2);
                    %     err(jj) = (err(jj) + ((vm./(vc.*Na)).*(z(kk) - qfun))^2);
                    % end
                end
            end
        end
    case 'STATZ2'
        % Discretize the data range into 'nbins' evenly spaced bins of pressure
        % ranges ('nbins = 1' for MLE')
        [bins] = discretize(x,nbins);
        % Create vector for storing sum of error for each bin
        err = zeros(nbins,1);
        Na = 6.022e20; % Avogadros constant [molecules/mmol]
        % Number of data points
        Nt = length(bins);
        omega = round(varargin{1}).*isoRef(1);
        beta = varargin{2}.*isoRef(2);
        b01NP = varargin{3}.*isoRef(3);
        b01LP = varargin{4}.*isoRef(3);
        delU1NP = varargin{5}.*isoRef(4);
        delU1LP = varargin{6}.*isoRef(4);
        delvc = varargin{7};
        vc = varargin{8};
        vm = varargin{9};
        betaVdW = varargin{10};
        vc_og = varargin{11};
        cutoffq = varargin{12};
        cutoffq1 = cutoffq(1);
        cutoffq2 = cutoffq(2);
        cutoffq0 = cutoffq(3);

        expData = [x,z,y];
        expData = sortrows(expData,2);
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

        xNP = [];
        yNP = [];
        zNP = [];

        xLP = [];
        yLP = [];
        zLP = [];
        if cutoffq2 == 55
            if cutoffq0.*((vc.*Na)./(vm)) > 0
                for ii = 1:length(temperatureValues)
                    indexCutoff2(ii,1) = qRefIndexTemp(ii,1);
                    if isempty(find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) < cutoffq0.*((vc.*Na)./(vm)),1,"last"))
                        indexCutoff2(ii,2) = qRefIndexTemp(ii,end-1);
                    else
                        indexCutoff2(ii,2) = indexCutoff2(ii,1) + find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) < cutoffq0.*((vc.*Na)./(vm)),1,"last");

                    end

                    if isempty(find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) < cutoffq1.*((vc.*Na)./(vm)),1,"last"))
                        xNP = [xNP; x(indexCutoff2(ii,2):indexCutoff2(ii,3))];
                        yNP = [yNP; y(indexCutoff2(ii,2):indexCutoff2(ii,3))];
                        zNP = [zNP; z(indexCutoff2(ii,2):indexCutoff2(ii,3))];
                        xLP = [xLP; x(indexCutoff2(ii,1):indexCutoff2(ii,2))];
                        yLP = [yLP; y(indexCutoff2(ii,1):indexCutoff2(ii,2))];
                        zLP = [zLP; z(indexCutoff2(ii,1):indexCutoff2(ii,2))];
                    else
                        indexCutoff2(ii,3) = indexCutoff2(ii,1) -1 + find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) < cutoffq1.*((vc.*Na)./(vm)),1,"last");
                        xNP = [xNP; x(indexCutoff2(ii,2):indexCutoff2(ii,3))];
                        yNP = [yNP; y(indexCutoff2(ii,2):indexCutoff2(ii,3))];
                        zNP = [zNP; z(indexCutoff2(ii,2):indexCutoff2(ii,3))];
                        xLP = [xLP; x(indexCutoff2(ii,3)+2:qRefIndexTemp(ii,2)); x(indexCutoff2(ii,1):indexCutoff2(ii,2))];
                        yLP = [yLP; y(indexCutoff2(ii,3)+2:qRefIndexTemp(ii,2)); y(indexCutoff2(ii,1):indexCutoff2(ii,2))];
                        zLP = [zLP; z(indexCutoff2(ii,3)+2:qRefIndexTemp(ii,2)); z(indexCutoff2(ii,1):indexCutoff2(ii,2))];
                    end
                end
            end
        else
            if cutoffq0.*((vc.*Na)./(vm)) > 0
                for ii = 1:length(temperatureValues)
                    indexCutoff2(ii,1) = qRefIndexTemp(ii,1);
                    if isempty(find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) < cutoffq0.*((vc.*Na)./(vm)),1,"last"))
                    else
                        indexCutoff2(ii,2) = indexCutoff2(ii,1) - 1 + find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) < cutoffq0.*((vc.*Na)./(vm)),1,"last");
                        xNP = [xNP; x(indexCutoff2(ii,1):indexCutoff2(ii,2))];
                        yNP = [yNP; y(indexCutoff2(ii,1):indexCutoff2(ii,2))];
                        zNP = [zNP; z(indexCutoff2(ii,1):indexCutoff2(ii,2))];
                    end

                    if isempty(find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) < cutoffq1.*((vc.*Na)./(vm)),1,"last"))
                    else
                        indexCutoff2(ii,3) = indexCutoff2(ii,1) + find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) < cutoffq1.*((vc.*Na)./(vm)),1,"last");
                        xLP = [xLP; x(indexCutoff2(ii,3):qRefIndexTemp(ii,2))];
                        yLP = [yLP; y(indexCutoff2(ii,3):qRefIndexTemp(ii,2))];
                        zLP = [zLP; z(indexCutoff2(ii,3):qRefIndexTemp(ii,2))];
                    end
                end
            end
        end
        % % Find qref for the experimental data
        % qRefIndexTempNP = [];
        % qRefIndexTempLP = [];
        % for ii = 1:length(temperatureValues)
        %     if ~isempty(find(yNP == temperatureValues(ii),1,'first')) &&
        %     qRefIndexTempNP(ii,1) = find(yNP == temperatureValues(ii),1,'first');
        %     qRefIndexTempNP(ii,2) = find(yNP == temperatureValues(ii),1,'last');
        %     end
        %     if ~isempty(find(yLP == temperatureValues(ii),1,'first'))
        %     qRefIndexTempLP(ii,1) = find(yLP == temperatureValues(ii),1,'first');
        %     qRefIndexTempLP(ii,2) = find(yLP == temperatureValues(ii),1,'last');
        %     end
        % end
        % 
        % qRefMaxNP = max(zNP(qRefIndexTempNP(:,2)));
        % qRefTempNP = zNP(qRefIndexTempNP(:,2));
        % normalizationFactorTempNP = qRefMaxNP./qRefTempNP;
        % if ~isempty(find(yLP == temperatureValues(ii),1,'first'))
        %     qRefMaxLP = max(zLP(qRefIndexTempLP(:,2)));
        %     qRefTempLP = zLP(qRefIndexTempLP(:,2));
        %     normalizationFactorTempLP = qRefMaxLP./qRefTempLP;
        %     normalizationFactorLP = zeros(length(xLP),1);
        % else
        %     normalizationFactorTempLP = ones(length(temperatureValues),1);
        %     normalizationFactorLP = zeros(length(xLP),1);
        % end
        % normalizationFactorNP = zeros(length(xNP),1);
        % 
        % for ii = 1:length(temperatureValues)
        %     normalizationFactorNP(qRefIndexTempNP(ii,1):qRefIndexTempNP(ii,2),1) = normalizationFactorTempNP(ii);
        %     if ~isempty(find(yLP == temperatureValues(ii),1,'first'))
        %         normalizationFactorLP(qRefIndexTempLP(ii,1):qRefIndexTempLP(ii,2),1) = normalizationFactorTempLP(ii);
        %     else
        % 
        %     end
        %     % normalizationFactor(qRefIndexTemp(ii,1):qRefIndexTemp(ii,1)+5,1) = normalizationFactorTemp(ii);
        % end

        % figure
        % scatter(xLP,zLP.*1,'ob')
        % hold on
        % scatter(xNP,zNP.*1,'or')
        % scatter(x,z.*1,'xk')
    case 'STATZ2gamma'
        % Discretize the data range into 'nbins' evenly spaced bins of pressure
        % ranges ('nbins = 1' for MLE')
        [bins] = discretize(x,nbins);
        % Create vector for storing sum of error for each bin
        err = zeros(nbins,1);
        Na = 6.022e20; % Avogadros constant [molecules/mmol]
        % Number of data points
        Nt = length(bins);
        omega = round(varargin{1}).*isoRef(1);
        beta = varargin{2}.*isoRef(2);
        b01NP = varargin{3}.*isoRef(3);
        b01LP = varargin{4}.*isoRef(3);
        delU1NP = varargin{5}.*isoRef(4);
        delU1LP = varargin{6}.*isoRef(4);
        delvc = varargin{7};
        gamma = varargin{8};
        vc = varargin{9};
        vm = varargin{10};
        betaVdW = varargin{11};
        vc_og = varargin{12};
        cutoffq = varargin{13};
        cutoffq1 = cutoffq(1);
        cutoffq2 = cutoffq(2);
        cutoffq0 = cutoffq(3);

        expData = [x,z,y];
        expData = sortrows(expData,2);
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

        xNP = [];
        yNP = [];
        zNP = [];

        xLP = [];
        yLP = [];
        zLP = [];
        if cutoffq2 == 55
            if cutoffq0.*((vc.*Na)./(vm)) > 0
                for ii = 1:length(temperatureValues)
                    indexCutoff2(ii,1) = qRefIndexTemp(ii,1);
                    if isempty(find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) < cutoffq0.*((vc.*Na)./(vm)),1,"last"))
                        indexCutoff2(ii,2) = qRefIndexTemp(ii,end-1);
                    else
                        indexCutoff2(ii,2) = indexCutoff2(ii,1) + find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) < cutoffq0.*((vc.*Na)./(vm)),1,"last");

                    end

                    if isempty(find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) < cutoffq1.*((vc.*Na)./(vm)),1,"last"))
                        xNP = [xNP; x(indexCutoff2(ii,2):indexCutoff2(ii,3))];
                        yNP = [yNP; y(indexCutoff2(ii,2):indexCutoff2(ii,3))];
                        zNP = [zNP; z(indexCutoff2(ii,2):indexCutoff2(ii,3))];
                        xLP = [xLP; x(indexCutoff2(ii,1):indexCutoff2(ii,2))];
                        yLP = [yLP; y(indexCutoff2(ii,1):indexCutoff2(ii,2))];
                        zLP = [zLP; z(indexCutoff2(ii,1):indexCutoff2(ii,2))];
                    else
                        indexCutoff2(ii,3) = indexCutoff2(ii,1) -1 + find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) < cutoffq1.*((vc.*Na)./(vm)),1,"last");
                        xNP = [xNP; x(indexCutoff2(ii,2):indexCutoff2(ii,3))];
                        yNP = [yNP; y(indexCutoff2(ii,2):indexCutoff2(ii,3))];
                        zNP = [zNP; z(indexCutoff2(ii,2):indexCutoff2(ii,3))];
                        xLP = [xLP; x(indexCutoff2(ii,3)+2:qRefIndexTemp(ii,2)); x(indexCutoff2(ii,1):indexCutoff2(ii,2))];
                        yLP = [yLP; y(indexCutoff2(ii,3)+2:qRefIndexTemp(ii,2)); y(indexCutoff2(ii,1):indexCutoff2(ii,2))];
                        zLP = [zLP; z(indexCutoff2(ii,3)+2:qRefIndexTemp(ii,2)); z(indexCutoff2(ii,1):indexCutoff2(ii,2))];
                    end
                end
            end
        else
            if cutoffq0.*((vc.*Na)./(vm)) > 0
                for ii = 1:length(temperatureValues)
                    indexCutoff2(ii,1) = qRefIndexTemp(ii,1);
                    if isempty(find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) < cutoffq0.*((vc.*Na)./(vm)),1,"last"))
                    else
                        indexCutoff2(ii,2) = indexCutoff2(ii,1) - 1 + find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) < cutoffq0.*((vc.*Na)./(vm)),1,"last");
                        xNP = [xNP; x(indexCutoff2(ii,1):indexCutoff2(ii,2))];
                        yNP = [yNP; y(indexCutoff2(ii,1):indexCutoff2(ii,2))];
                        zNP = [zNP; z(indexCutoff2(ii,1):indexCutoff2(ii,2))];
                    end

                    if isempty(find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) < cutoffq1.*((vc.*Na)./(vm)),1,"last"))
                    else
                        indexCutoff2(ii,3) = indexCutoff2(ii,1) + find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) < cutoffq1.*((vc.*Na)./(vm)),1,"last");
                        xLP = [xLP; x(indexCutoff2(ii,3):qRefIndexTemp(ii,2))];
                        yLP = [yLP; y(indexCutoff2(ii,3):qRefIndexTemp(ii,2))];
                        zLP = [zLP; z(indexCutoff2(ii,3):qRefIndexTemp(ii,2))];
                    end
                end
            end
        end
        % % Find qref for the experimental data
        % qRefIndexTempNP = [];
        % qRefIndexTempLP = [];
        % for ii = 1:length(temperatureValues)
        %     if ~isempty(find(yNP == temperatureValues(ii),1,'first')) &&
        %     qRefIndexTempNP(ii,1) = find(yNP == temperatureValues(ii),1,'first');
        %     qRefIndexTempNP(ii,2) = find(yNP == temperatureValues(ii),1,'last');
        %     end
        %     if ~isempty(find(yLP == temperatureValues(ii),1,'first'))
        %     qRefIndexTempLP(ii,1) = find(yLP == temperatureValues(ii),1,'first');
        %     qRefIndexTempLP(ii,2) = find(yLP == temperatureValues(ii),1,'last');
        %     end
        % end
        % 
        % qRefMaxNP = max(zNP(qRefIndexTempNP(:,2)));
        % qRefTempNP = zNP(qRefIndexTempNP(:,2));
        % normalizationFactorTempNP = qRefMaxNP./qRefTempNP;
        % if ~isempty(find(yLP == temperatureValues(ii),1,'first'))
        %     qRefMaxLP = max(zLP(qRefIndexTempLP(:,2)));
        %     qRefTempLP = zLP(qRefIndexTempLP(:,2));
        %     normalizationFactorTempLP = qRefMaxLP./qRefTempLP;
        %     normalizationFactorLP = zeros(length(xLP),1);
        % else
        %     normalizationFactorTempLP = ones(length(temperatureValues),1);
        %     normalizationFactorLP = zeros(length(xLP),1);
        % end
        % normalizationFactorNP = zeros(length(xNP),1);
        % 
        % for ii = 1:length(temperatureValues)
        %     normalizationFactorNP(qRefIndexTempNP(ii,1):qRefIndexTempNP(ii,2),1) = normalizationFactorTempNP(ii);
        %     if ~isempty(find(yLP == temperatureValues(ii),1,'first'))
        %         normalizationFactorLP(qRefIndexTempLP(ii,1):qRefIndexTempLP(ii,2),1) = normalizationFactorTempLP(ii);
        %     else
        % 
        %     end
        %     % normalizationFactor(qRefIndexTemp(ii,1):qRefIndexTemp(ii,1)+5,1) = normalizationFactorTemp(ii);
        % end

        % figure
        % scatter(xLP,zLP.*1,'ob')
        % hold on
        % scatter(xNP,zNP.*1,'or')
        % scatter(x,z.*1,'xk')
        
        for jj = 1:length(err)
            for kk = 1:length(zNP)
                % if omega > vc/beta
                %     qfun = 0;
                % else
                qfunNP = computeStatSTALoading2(xNP(kk),yNP(kk),b01NP,b01NP,delU1NP,delU1NP,beta,0,0,0,floor((vc.*(1-delvc))./betaVdW),vc.*(1-delvc),gamma);
                
                % end
                if bins(kk) == jj
                    % if kk < round(length(err)./10)
                    % err(jj) = (err(jj) + (normalizationFactor(kk).*(vm./(vc.*Na)).*(z(kk) - qfun))^2);
                    % err(jj) = (err(jj) + ((zNP(kk)./((vc.*Na)./(vm)) - qfunNP./((vc.*(1-delvc).*Na)./(vm))))^2);
                    if isnan(qfunNP)
                        err(jj) = 1000000;
                    else
                        % err(jj) = (err(jj) + ((zNP(kk)./((vc_og.*Na)./(vm)) - qfunNP./((vc.*(1-delvc).*Na)./(vm))))^2);
                        err(jj) = (err(jj) + ((zNP(kk) - qfunNP))^2);
                    end
                    % else
                    % err(jj) = (err(jj) + (normalizationFactor(kk).*(vm./(vc.*Na)).*(z(kk) - qfun))^2);
                    %     err(jj) = (err(jj) + ((vm./(vc.*Na)).*(z(kk) - qfun))^2);
                    % end
                end
            end
            for kk = 1:length(zLP)
                % if omega > vc/beta
                %     qfun = 0;
                % else
                qfunLP = computeStatSTALoading2(xLP(kk),yLP(kk),b01LP,b01LP,delU1LP,delU1LP,beta,0,0,0,omega,vc,gamma);
                % end
                if bins(kk) == jj
                    % if kk < round(length(err)./10)
                    % err(jj) = (err(jj) + (normalizationFactor(kk).*(vm./(vc.*Na)).*(z(kk) - qfun))^2);
                    % err(jj) = (err(jj) + ((zLP(kk)./((vc.*Na)./(vm)) - qfunLP./((vc.*Na)./(vm))))^2);
                    % err(jj) = (err(jj) + ((zNP(kk)./((vc.*Na)./(vm)) - qfunNP./((vc.*(1-delvc).*Na)./(vm))))^2);
                    if isnan(qfunLP)
                        err(jj) = 1000000;
                    else                 
                        % err(jj) = (err(jj) + ((zLP(kk)./((vc_og.*Na)./(vm)) - qfunLP./((vc.*Na)./(vm))))^2);
                        err(jj) = (err(jj) + ((zLP(kk) - qfunLP))^2);
                    end
                    % else
                    % err(jj) = (err(jj) + (normalizationFactor(kk).*(vm./(vc.*Na)).*(z(kk) - qfun))^2);
                    %     err(jj) = (err(jj) + ((vm./(vc.*Na)).*(z(kk) - qfun))^2);
                    % end
                end
            end
        end
        1+1;
    case 'SSL2'
        % Discretize the data range into 'nbins' evenly spaced bins of pressure
        % ranges ('nbins = 1' for MLE')
        [bins] = discretize(x,nbins);
        % Create vector for storing sum of error for each bin
        err = zeros(nbins,1);
        Na = 6.022e20; % Avogadros constant [molecules/mmol]
        % Number of data points
        Nt = length(bins);
        qsNP = varargin{1}.*isoRef(1);
        b01NP = varargin{2}.*isoRef(2);
        delU1NP = varargin{3}.*isoRef(3);
        qsLP = varargin{4}.*isoRef(4);
        b01LP = varargin{5}.*isoRef(5);
        delU1LP = varargin{6}.*isoRef(6);
        cutoffq = varargin{7};
        cutoffq1 = cutoffq(1);
        cutoffq2 = cutoffq(2);
        cutoffq0 = cutoffq(3);

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

        xNP = [];
        yNP = [];
        zNP = [];

        xLP = [];
        yLP = [];
        zLP = [];

        if cutoffq0 > 0
            for ii = 1:length(temperatureValues)
                indexCutoff2(ii,1) = qRefIndexTemp(ii,1);
                indexCutoff2(ii,2) = indexCutoff2(ii,1) + find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) < cutoffq0,1,"last");
                xNP = [xNP; x(indexCutoff2(ii,1):indexCutoff2(ii,2))];
                yNP = [yNP; y(indexCutoff2(ii,1):indexCutoff2(ii,2))];
                zNP = [zNP; z(indexCutoff2(ii,1):indexCutoff2(ii,2))];

                if isempty(find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) < cutoffq1,1,"last"))
                else
                    indexCutoff2(ii,3) = indexCutoff2(ii,1) + find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) < cutoffq1,1,"last");
                    xLP = [xLP; x(indexCutoff2(ii,3):qRefIndexTemp(ii,2))];
                    yLP = [yLP; y(indexCutoff2(ii,3):qRefIndexTemp(ii,2))];
                    zLP = [zLP; z(indexCutoff2(ii,3):qRefIndexTemp(ii,2))];
                end
            end
        else

        end

        % Find qref for the experimental data
        qRefIndexTempNP = zeros(length(temperatureValues),1);
        qRefIndexTempLP = zeros(length(temperatureValues),1);
        for ii = 1:length(temperatureValues)
            qRefIndexTempNP(ii,1) = find(yNP == temperatureValues(ii),1,'first');
            qRefIndexTempNP(ii,2) = find(yNP == temperatureValues(ii),1,'last');
            if ~isempty(find(yLP == temperatureValues(ii),1,'first'))
            qRefIndexTempLP(ii,1) = find(yLP == temperatureValues(ii),1,'first');
            qRefIndexTempLP(ii,2) = find(yLP == temperatureValues(ii),1,'last');
            else

            end
        end

        qRefMaxNP = max(zNP(qRefIndexTempNP(:,2)));
        qRefTempNP = zNP(qRefIndexTempNP(:,2));
        normalizationFactorTempNP = qRefMaxNP./qRefTempNP;
        if ~isempty(find(yLP == temperatureValues(ii),1,'first'))
            qRefMaxLP = max(zLP(qRefIndexTempLP(:,2)));
            qRefTempLP = zLP(qRefIndexTempLP(:,2));
            normalizationFactorTempLP = qRefMaxLP./qRefTempLP;
            normalizationFactorLP = zeros(length(xLP),1);
        else
            normalizationFactorTempLP = ones(length(temperatureValues),1);
            normalizationFactorLP = zeros(length(xLP),1);
        end
        normalizationFactorNP = zeros(length(xNP),1);

        for ii = 1:length(temperatureValues)
            normalizationFactorNP(qRefIndexTempNP(ii,1):qRefIndexTempNP(ii,2),1) = normalizationFactorTempNP(ii);
            if ~isempty(find(yLP == temperatureValues(ii),1,'first'))
                normalizationFactorLP(qRefIndexTempLP(ii,1):qRefIndexTempLP(ii,2),1) = normalizationFactorTempLP(ii);
            else

            end
            % normalizationFactor(qRefIndexTemp(ii,1):qRefIndexTemp(ii,1)+5,1) = normalizationFactorTemp(ii);
        end

        % figure
        % scatter(xLP,zLP.*1)
        % hold on
        % scatter(xNP,zNP.*1)

        
        for jj = 1:length(err)
            for kk = 1:length(zNP)
                % if omega > vc/beta
                %     qfun = 0;
                % else
                qfunNP = qsNP*(b01NP*xNP(kk)*exp(delU1NP/(8.314*yNP(kk))))/(1+(b01NP*xNP(kk)*exp(delU1NP/(8.314*yNP(kk)))));
                % end
                if bins(kk) == jj
                    if isnan(qfunNP)
                        err(jj) = 1000000;
                    else
                        err(jj) = (err(jj) + ((zNP(kk) - qfunNP))^2);
                    end

                end
            end
            for kk = 1:length(zLP)
                qfunLP = qsLP*(b01LP*xLP(kk)*exp(delU1LP/(8.314*yLP(kk))))/(1+(b01LP*xLP(kk)*exp(delU1LP/(8.314*yLP(kk)))));
                if bins(kk) == jj
                    if isnan(qfunLP)
                        err(jj) = 1000000;
                    else                 
                        err(jj) = (err(jj) + ((zLP(kk) - qfunLP))^2);
                    end
                end
            end
        end
        1+1;        
    case 'STATZE'
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
        ebyk = varargin{5}.*isoRef(5);
        vc = varargin{6};
        vm = varargin{7};
        Na = 6.022e20; % Avogadros constant [molecules/mmol]
        for jj = 1:length(err)
            for kk = 1:length(bins)
                % if omega > vc/beta
                %     qfun = 0;
                % else
                    qfun = computeStatZELoading(x(kk),y(kk),b01,delU1,beta,omega,ebyk,vc);
                % end
                if bins(kk) == jj
                    err(jj) = (err(jj) + (normalizationFactor(kk).*(vm./(vc.*Na)).*(z(kk) - qfun))^2);
                end
            end
        end
    case 'STATZGO'
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
        beta =  varargin{2}.*isoRef(2);
        b01 =   varargin{3}.*isoRef(3);
        delU1 = varargin{4}.*isoRef(4);
        delU2 = varargin{5}.*isoRef(5);
        kgate = varargin{6}.*isoRef(6);
        cgate = varargin{7}.*isoRef(7);
        vc = varargin{8};
        vm = varargin{9};
        Na = 6.022e20; % Avogadros constant [molecules/mmol]
        for jj = 1:length(err)
            for kk = 1:length(bins)
                % if omega > vc/beta
                %     qfun = 0;
                % else
                    qfun = computeStatZGATELoading2(x(kk),y(kk),b01,delU1,delU2,beta,kgate,cgate,omega,vc);
                % end
                if bins(kk) == jj
                    err(jj) = (err(jj) + ((z(kk) - qfun))^2);
                end
            end
        end
    case 'SSLSTA'
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

        qsNP = varargin{1}.*isoRef(1);
        b01NP = varargin{2}.*isoRef(2);
        delU1NP = varargin{3}.*isoRef(3);
        qsLP = varargin{4}.*isoRef(4);
        b01LP = varargin{5}.*isoRef(5);
        delU1LP = varargin{6}.*isoRef(6);
        kgate = varargin{7}.*isoRef(7);
        cgate = varargin{8}.*isoRef(8);
        sval =  varargin{9}.*isoRef(9);

        for jj = 1:length(err)
            for kk = 1:length(bins)
                % if omega > vc/beta
                %     qfun = 0;
                % else
                yval = ((1+(b01LP*x(kk)*exp(delU1LP/(8.314*y(kk))))).^qsLP)./((1+(b01NP*x(kk)*exp(delU1NP/(8.314*y(kk))))).^qsNP).* ...
                    exp(-(kgate-y(kk).*cgate)./(8.314.*y(kk)));
                sigmaval = yval.^sval./(1+yval.^sval);
                qfun = (1-sigmaval).*qsNP*(b01NP*x(kk)*exp(delU1NP/(8.314*y(kk))))/(1+(b01NP*x(kk)*exp(delU1NP/(8.314*y(kk))))) + ...
                    sigmaval.*qsLP*(b01LP*x(kk)*exp(delU1LP/(8.314*y(kk))))/(1+(b01LP*x(kk)*exp(delU1LP/(8.314*y(kk)))));
                % end
                if bins(kk) == jj
                    if isnan(qfun)
                        err(jj) = err(jj) + 1000;
                    else
                        err(jj) = normalizationFactor(kk).*(err(jj) + ((z(kk) - qfun))^2);
                    end
                end
            end
        end
    case 'STATSTA2'
        % Discretize the data range into 'nbins' evenly spaced bins of pressure
        % ranges ('nbins = 1' for MLE')
        [bins] = discretize(x,nbins);
        % Create vector for storing sum of error for each bin
        err = zeros(nbins,1);
        % Number of data points
        Nt = length(bins);

        expData = [x,z,y];
        expData = sortrows(expData,3);

        

        omega = round(varargin{1}).*isoRef(1);
        beta =  varargin{2}.*isoRef(2);
        b01 =   varargin{3}.*isoRef(3);
        b02 =   varargin{4}.*isoRef(3);
        delU1 = varargin{5}.*isoRef(4);
        delU2 = varargin{6}.*isoRef(4);
        kgate = varargin{7}.*isoRef(6);
        cgate = varargin{8}.*isoRef(7);
        delvc = varargin{9}.*isoRef(8);
        vc = varargin{10};
        vm = varargin{11};
        vc_og = varargin{12};
        Na = 6.022e20; % Avogadros constant [molecules/mmol]
        cutoffq = varargin{13};
        cutoffq1 = cutoffq(1).*((vc_og.*Na)./(vm));
        cutoffq2 = cutoffq(2).*((vc_og.*Na)./(vm));
        cutoffq0 = cutoffq(3).*((vc_og.*Na)./(vm));

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
        
        indexCutoff = zeros(length(temperatureValues),2);

        % indexCutoff2 = zeros(length(temperatureValues),2);
        % x_new1 = [];
        % y_new1 = [];
        % z_new1 = [];
        % if cutoffq2 > 0
        %     for ii = 1:length(temperatureValues)
        %         indexCutoff2(ii,1) = qRefIndexTemp(ii,1);
        %         indexCutoff2(ii,2) = indexCutoff2(ii,1) + find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) > cutoffq2,1,"first");
        %         x_new1 = [x_new1; x(indexCutoff2(ii,1):indexCutoff2(ii,2))];
        %         y_new1 = [y_new1; y(indexCutoff2(ii,1):indexCutoff2(ii,2))];
        %         z_new1 = [z_new1; z(indexCutoff2(ii,1):indexCutoff2(ii,2))];
        %         normalizationFactor(indexCutoff(ii,2)+1:qRefIndexTemp(ii,2),1) = normalizationFactor(indexCutoff(ii,2)+1:qRefIndexTemp(ii,2),1)./max(z(indexCutoff(ii,1):indexCutoff(ii,2)));
        %     end
        % else
            x_new1 = x;
            y_new1 = y;
            z_new1 = z;
        % end

        % temperatureValues = unique(y_new1);
        % qRefIndexTemp = zeros(length(temperatureValues),1);
        % for ii = 1:length(temperatureValues)
        %     qRefIndexTemp(ii,1) = find(y_new1 == temperatureValues(ii),1,'first');
        %     qRefIndexTemp(ii,2) = find(y_new1 == temperatureValues(ii),1,'last');
        % end

        % indexCutoff3 = zeros(length(temperatureValues),2);
        x_new = x_new1;
        y_new = y_new1;
        z_new = z_new1;

        % if cutoffq0 == 0
            for ii = 1:length(temperatureValues)
                % indexCutoff3(ii,2) = qRefIndexTemp(ii,2);
                % indexCutoff3(ii,1) = indexCutoff3(ii,2) - find(z_new1(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) > cutoffq0,1,"first");
                % x_new = [x_new; x_new1(indexCutoff3(ii,1):indexCutoff3(ii,2))];
                % y_new = [y_new; y_new1(indexCutoff3(ii,1):indexCutoff2(ii,2))];
                % z_new = [z_new; z_new1(indexCutoff3(ii,1):indexCutoff2(ii,2))];
                normalizationFactor(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) = normalizationFactorTemp(ii);
            end
        % else
            % normalizationFactor(:) = 1;
        % end

        % temperatureValues = unique(y_new);
        % qRefIndexTemp = zeros(length(temperatureValues),1);
        % for ii = 1:length(temperatureValues)
        %     qRefIndexTemp(ii,1) = find(y_new == temperatureValues(ii),1,'first');
        %     qRefIndexTemp(ii,2) = find(y_new == temperatureValues(ii),1,'last');
        % end

        
        
        % normalizationFactor = zeros(length(z_new),1);
        % normalizationFactor(:) = 1;

        if cutoffq0 > 0
            for ii = 1:length(temperatureValues)
                indexCutoff(ii,1) = qRefIndexTemp(ii,1);
                if ~isempty(find(z_new(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) > cutoffq0,1,"first"))
                    indexCutoff(ii,2) = indexCutoff(ii,1) -1 + find(z_new(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) > cutoffq0,1,"first");
                else
                    indexCutoff(ii,2) = qRefIndexTemp(ii,2);
                end
                normalizationFactor(indexCutoff(ii,1):indexCutoff(ii,2)) = 1.*normalizationFactor(indexCutoff(ii,1):indexCutoff(ii,2));
                normalizationFactor(indexCutoff(ii,2):qRefIndexTemp(ii,2)) = 1.*normalizationFactor(indexCutoff(ii,2):qRefIndexTemp(ii,2));
                % normalizationFactor(indexCutoff(ii,1):indexCutoff(ii,2)) = 1.*normalizationFactor(indexCutoff(ii,1):indexCutoff(ii,2))./length(normalizationFactor(indexCutoff(ii,1):indexCutoff(ii,2)));
                % normalizationFactor(indexCutoff(ii,2):qRefIndexTemp(ii,2)) = 1.*normalizationFactor(indexCutoff(ii,2):qRefIndexTemp(ii,2))./length(normalizationFactor(indexCutoff(ii,2):qRefIndexTemp(ii,2)));
            end
        else
            normalizationFactor(:) = 1;
        end


        % figure
        % scatter(x_new,z_new)
        % hold on
        % yline(cutoffq0)
        % set(gca,'XScale','log')

        [bins] = discretize(x_new,nbins);
        % Create vector for storing sum of error for each bin
        err = zeros(nbins,1);
        % Number of data points
        Nt = length(bins);


        for jj = 1:length(err)
            for kk = 1:length(bins)
                % if omega > vc/beta
                %     qfun = 0;
                % else
                qfun = computeStatSTALoading2(x_new(kk),y_new(kk),b01,b02,delU1,delU2,beta,kgate,cgate,delvc,omega,vc);
                % end
                if bins(kk) == jj
                    if isnan(qfun)
                        err(jj) = err(jj) + 1e50;
                        % break
                    else
                        % if cutoffq2./((vc_og.*Na)./(vm)) == 99 && kgate-cgate.*y_new(kk) > 0
                        %     err(jj) = err(jj) + 1e50;
                        % elseif round(cutoffq2./((vc_og.*Na)./(vm))) == 98 && cgate < 0
                        %     err(jj) = err(jj) + 1e50;
                        %     % break
                        % else
                            % err(jj) = (err(jj) + (normalizationFactor(kk).*(z_new(kk)./((vc_og.*Na)./(vm)) - qfun./((vc.*Na)./(vm))))^2);
                            % err(jj) = (err(jj) + normalizationFactor(kk).*((z_new(kk)./((vc_og.*Na)./(vm)) - qfun./((vc.*Na)./(vm))))^2);
                            if cutoffq0 > 0
                                err(jj) = (err(jj) + (normalizationFactor(kk).*(z_new(kk) - qfun))^2);
                            else
                                err(jj) = (err(jj) + (normalizationFactor(kk).*(z_new(kk) - qfun)./((vc_og.*Na)./(vm)))^2);
                            end
                        % end
                    end

                end
            end
        end
    case 'STATSTA2gamma'
        % Discretize the data range into 'nbins' evenly spaced bins of pressure
        % ranges ('nbins = 1' for MLE')
        [bins] = discretize(x,nbins);
        % Create vector for storing sum of error for each bin
        err = zeros(nbins,1);
        % Number of data points
        Nt = length(bins);

        expData = [x,z,y];
        expData = sortrows(expData,3);

        

        omega = round(varargin{1}).*isoRef(1);
        beta =  varargin{2}.*isoRef(2);
        b01 =   varargin{3}.*isoRef(3);
        b02 =   varargin{4}.*isoRef(3);
        delU1 = varargin{5}.*isoRef(4);
        delU2 = varargin{6}.*isoRef(4);
        kgate = varargin{7}.*isoRef(6);
        cgate = varargin{8}.*isoRef(7);
        delvc = varargin{9}.*isoRef(8);
        gamma = varargin{10}.*isoRef(9);
        vc = varargin{11};
        vm = varargin{12};
        vc_og = varargin{13};
        Na = 6.022e20; % Avogadros constant [molecules/mmol]
        cutoffq = varargin{14};
        cutoffq1 = cutoffq(1).*((vc_og.*Na)./(vm));
        cutoffq2 = cutoffq(2).*((vc_og.*Na)./(vm));
        cutoffq0 = cutoffq(3).*((vc_og.*Na)./(vm));

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
        
        indexCutoff = zeros(length(temperatureValues),2);

        % indexCutoff2 = zeros(length(temperatureValues),2);
        % x_new1 = [];
        % y_new1 = [];
        % z_new1 = [];
        % if cutoffq2 > 0
        %     for ii = 1:length(temperatureValues)
        %         indexCutoff2(ii,1) = qRefIndexTemp(ii,1);
        %         indexCutoff2(ii,2) = indexCutoff2(ii,1) + find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) > cutoffq2,1,"first");
        %         x_new1 = [x_new1; x(indexCutoff2(ii,1):indexCutoff2(ii,2))];
        %         y_new1 = [y_new1; y(indexCutoff2(ii,1):indexCutoff2(ii,2))];
        %         z_new1 = [z_new1; z(indexCutoff2(ii,1):indexCutoff2(ii,2))];
        %         normalizationFactor(indexCutoff(ii,2)+1:qRefIndexTemp(ii,2),1) = normalizationFactor(indexCutoff(ii,2)+1:qRefIndexTemp(ii,2),1)./max(z(indexCutoff(ii,1):indexCutoff(ii,2)));
        %     end
        % else
            x_new1 = x;
            y_new1 = y;
            z_new1 = z;
        % end

        % temperatureValues = unique(y_new1);
        % qRefIndexTemp = zeros(length(temperatureValues),1);
        % for ii = 1:length(temperatureValues)
        %     qRefIndexTemp(ii,1) = find(y_new1 == temperatureValues(ii),1,'first');
        %     qRefIndexTemp(ii,2) = find(y_new1 == temperatureValues(ii),1,'last');
        % end

        % indexCutoff3 = zeros(length(temperatureValues),2);
        x_new = x_new1;
        y_new = y_new1;
        z_new = z_new1;

        % if cutoffq0 == 0
            for ii = 1:length(temperatureValues)
                % indexCutoff3(ii,2) = qRefIndexTemp(ii,2);
                % indexCutoff3(ii,1) = indexCutoff3(ii,2) - find(z_new1(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) > cutoffq0,1,"first");
                % x_new = [x_new; x_new1(indexCutoff3(ii,1):indexCutoff3(ii,2))];
                % y_new = [y_new; y_new1(indexCutoff3(ii,1):indexCutoff2(ii,2))];
                % z_new = [z_new; z_new1(indexCutoff3(ii,1):indexCutoff2(ii,2))];
                normalizationFactor(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) = normalizationFactorTemp(ii);
            end
        % else
            % normalizationFactor(:) = 1;
        % end

        % temperatureValues = unique(y_new);
        % qRefIndexTemp = zeros(length(temperatureValues),1);
        % for ii = 1:length(temperatureValues)
        %     qRefIndexTemp(ii,1) = find(y_new == temperatureValues(ii),1,'first');
        %     qRefIndexTemp(ii,2) = find(y_new == temperatureValues(ii),1,'last');
        % end

        
        
        % normalizationFactor = zeros(length(z_new),1);
        % normalizationFactor(:) = 1;

        if cutoffq0 > 0
            for ii = 1:length(temperatureValues)
                indexCutoff(ii,1) = qRefIndexTemp(ii,1);
                if ~isempty(find(z_new(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) > cutoffq0,1,"first"))
                    indexCutoff(ii,2) = indexCutoff(ii,1) -1 + find(z_new(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) > cutoffq0,1,"first");
                else
                    indexCutoff(ii,2) = qRefIndexTemp(ii,2);
                end
                normalizationFactor(indexCutoff(ii,1):indexCutoff(ii,2)) = 1.*normalizationFactor(indexCutoff(ii,1):indexCutoff(ii,2));
                normalizationFactor(indexCutoff(ii,2):qRefIndexTemp(ii,2)) = 1.*normalizationFactor(indexCutoff(ii,2):qRefIndexTemp(ii,2));
                % normalizationFactor(indexCutoff(ii,1):indexCutoff(ii,2)) = 1.*normalizationFactor(indexCutoff(ii,1):indexCutoff(ii,2))./length(normalizationFactor(indexCutoff(ii,1):indexCutoff(ii,2)));
                % normalizationFactor(indexCutoff(ii,2):qRefIndexTemp(ii,2)) = 1.*normalizationFactor(indexCutoff(ii,2):qRefIndexTemp(ii,2))./length(normalizationFactor(indexCutoff(ii,2):qRefIndexTemp(ii,2)));
            end
        else
            normalizationFactor(:) = 1;
        end


        % figure
        % scatter(x_new,z_new)
        % hold on
        % yline(cutoffq0)
        % set(gca,'XScale','log')

        [bins] = discretize(x_new,nbins);
        % Create vector for storing sum of error for each bin
        err = zeros(nbins,1);
        % Number of data points
        Nt = length(bins);


        for jj = 1:length(err)
            for kk = 1:length(bins)
                % if omega > vc/beta
                %     qfun = 0;
                % else
                qfun = computeStatSTALoading2(x_new(kk),y_new(kk),b01,b02,delU1,delU2,beta,kgate,cgate,delvc,omega,vc,gamma);
                % end
                if bins(kk) == jj
                    if isnan(qfun)
                        err(jj) = err(jj) + 1e50;
                        % break
                    else
                        % if cutoffq2./((vc_og.*Na)./(vm)) == 99 && kgate-cgate.*y_new(kk) > 0
                        %     err(jj) = err(jj) + 1e50;
                        % elseif round(cutoffq2./((vc_og.*Na)./(vm))) == 98 && cgate < 0
                        %     err(jj) = err(jj) + 1e50;
                        %     % break
                        % else
                            % err(jj) = (err(jj) + (normalizationFactor(kk).*(z_new(kk)./((vc_og.*Na)./(vm)) - qfun./((vc.*Na)./(vm))))^2);
                            % err(jj) = (err(jj) + normalizationFactor(kk).*((z_new(kk)./((vc_og.*Na)./(vm)) - qfun./((vc.*Na)./(vm))))^2);
                            if cutoffq0 > 0
                                err(jj) = (err(jj) + (normalizationFactor(kk).*(z_new(kk) - qfun))^2);
                            else
                                err(jj) = (err(jj) + (normalizationFactor(kk).*(z_new(kk) - qfun)./((vc_og.*Na)./(vm)))^2);
                            end
                        % end
                    end

                end
            end
        end
    case 'STATGO3'
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

        omega = round(varargin{1}).*isoRef(1);
        beta =  varargin{2}.*isoRef(2);
        b01 =   varargin{3}.*isoRef(3);
        delU1 = varargin{4}.*isoRef(4);
        delU2 = varargin{5}.*isoRef(5);
        kgate = varargin{6}.*isoRef(6);
        cgate = varargin{7}.*isoRef(7);
        delvc = varargin{8}.*isoRef(8);
        b02 = varargin{9}.*isoRef(9);
        vc = varargin{10};
        vm = varargin{11};
        cutoffq = varargin{12};

        Na = 6.022e20; % Avogadros constant [molecules/mmol]
        cutoffq1 = cutoffq(1).*((vc.*Na)./(vm));
        cutoffq2 = cutoffq(2).*((vc.*Na)./(vm));
        cutoffq0 = cutoffq(3).*((vc.*Na)./(vm));

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
        normalizationFactor = zeros(length(z),1);
        

        x_new = [];
        y_new = [];
        z_new = [];

        indexCutoff3 = zeros(length(temperatureValues),2);

        if cutoffq0 > 0
            for ii = 1:length(temperatureValues)
                normalizationFactor(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) = normalizationFactorTemp(ii);
                indexCutoff2(ii,1) = qRefIndexTemp(ii,1);
                indexCutoff2(ii,2) = indexCutoff2(ii,1) + find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) > cutoffq0,1,"first");
                % x_new1 = [x_new1; x(indexCutoff2(ii,1):indexCutoff2(ii,2))];
                % y_new1 = [y_new1; y(indexCutoff2(ii,1):indexCutoff2(ii,2))];
                % z_new1 = [z_new1; z(indexCutoff2(ii,1):indexCutoff2(ii,2))];
                normalizationFactor(indexCutoff2(ii,1):indexCutoff2(ii,2),1) = 1*normalizationFactor(indexCutoff2(ii,1):indexCutoff2(ii,2),1)./max(z(indexCutoff2(ii,1):indexCutoff2(ii,2)));
                normalizationFactor(indexCutoff2(ii,2)+1:qRefIndexTemp(ii,2),1) = 2*normalizationFactor(indexCutoff2(ii,2)+1:qRefIndexTemp(ii,2),1)./max(z(indexCutoff2(ii,2)+1:qRefIndexTemp(ii,2)));
            end
        else
            
        end

        x_new = x;
        y_new = y;
        z_new = z;

        % indexCutoff3 = zeros(length(temperatureValues),2);
        % x_new = [];
        % y_new = [];
        % z_new = [];
        % if cutoffq0 > 0
        %     for ii = 1:length(temperatureValues)
        %         indexCutoff3(ii,2) = qRefIndexTemp(ii,2);
        %         indexCutoff3(ii,1) = indexCutoff3(ii,2) - find(z_new1(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) > cutoffq0,1,"first");
        %         x_new = [x_new; x_new1(indexCutoff3(ii,1):indexCutoff3(ii,2))];
        %         y_new = [y_new; y_new1(indexCutoff3(ii,1):indexCutoff2(ii,2))];
        %         z_new = [z_new; z_new1(indexCutoff3(ii,1):indexCutoff2(ii,2))];
        %         % normalizationFactor(indexCutoff(ii,2)+1:qRefIndexTemp(ii,2),1) = normalizationFactor(indexCutoff(ii,2)+1:qRefIndexTemp(ii,2),1)./max(z(indexCutoff(ii,1):indexCutoff(ii,2)));
        %     end
        % else
        %     x_new = x_new1;
        %     y_new = y_new1;
        %     z_new = z_new1;
        % end
        % 
        % temperatureValues = unique(y_new);
        % qRefIndexTemp = zeros(length(temperatureValues),1);
        % for ii = 1:length(temperatureValues)
        %     qRefIndexTemp(ii,1) = find(y_new == temperatureValues(ii),1,'first');
        %     qRefIndexTemp(ii,2) = find(y_new == temperatureValues(ii),1,'last');
        % end

        % Find qref for the experimental data
        % qRefMax = max(z_new(qRefIndexTemp(:,2)));
        % qRefTemp = z_new(qRefIndexTemp(:,2));
        % normalizationFactorTemp = qRefMax./qRefTemp;
        % normalizationFactor = zeros(length(x_new),1);

        % if cutoffq1 > 0
        %     for ii = 1:length(temperatureValues)
        %         % indexCutoff(ii,1) = qRefIndexTemp(ii,1);
        %         % indexCutoff(ii,2) = indexCutoff(ii,1) + find(z_new(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) > cutoffq1,1,"first");
        %         % normalizationFactor(indexCutoff(ii,1):indexCutoff(ii,2),1)     = max(z_new(indexCutoff(ii,1):indexCutoff(ii,2),1))./z_new(indexCutoff(ii,1):indexCutoff(ii,2),1);
        %         % normalizationFactor(indexCutoff(ii,2)+1:qRefIndexTemp(ii,2),1) = ((z_new(indexCutoff(ii,2)+1:qRefIndexTemp(ii,2))))./max(z_new(indexCutoff(ii,2)+1:qRefIndexTemp(ii,2),1));
        %         % normalizationFactor(indexCutoff(ii,1):indexCutoff(ii,2),1) = normalizationFactor(indexCutoff(ii,1):indexCutoff(ii,2),1)./max(z_new(indexCutoff(ii,1):indexCutoff(ii,2)));
        %         % normalizationFactor(indexCutoff(ii,2)+1:qRefIndexTemp(ii,2),1) = normalizationFactor(indexCutoff(ii,2)+1:qRefIndexTemp(ii,2),1)./max(z(indexCutoff(ii,1):indexCutoff(ii,2)));
        %     end
        % else
        %     normalizationFactor(:) = 1;
        %     % normalizationFactor(qRefIndexTemp(:,2)-10:qRefIndexTemp(:,2)) = 4;
        % end

        [bins] = discretize(x_new,nbins);
        % Create vector for storing sum of error for each bin
        err = zeros(nbins,1);
        % Number of data points
        Nt = length(bins);

        % scatter(x, normalizationFactor.*z_new)

        
        for jj = 1:length(err)
            for kk = 1:length(bins)
                % if omega > vc/beta
                %     qfun = 0;
                % else
                    qfun = computeStatZGATELoadingGO(x_new(kk),y_new(kk),b01,b02,delU1,delU2,beta,kgate,cgate,delvc,omega,vc);
                % end
                if bins(kk) == jj
                    % if isnan(qfun)
                    %     err(jj) = err(jj) + 1000;
                    % else
                        err(jj) = (err(jj) + (normalizationFactor(kk).*(z_new(kk) - qfun))^2);
                    % end
                end
            end
        end
        % Calculate error for Statistical isotherm model for zeolites
    case 'STATZSips'
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
        gamma = varargin{5}.*isoRef(5);
        vc = varargin{6};
        vm = varargin{7};
        Na = 6.022e20; % Avogadros constant [molecules/mmol]
        for jj = 1:length(err)
            for kk = 1:length(bins)
                if omega > vc/beta
                    qfun = 0;
                else
                    qfun = computeStatZSipsLoading(x(kk),y(kk),b01,delU1,beta,omega,gamma,vc);
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