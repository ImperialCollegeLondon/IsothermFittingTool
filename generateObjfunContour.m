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
% Plot the contour plots for the objective function (MLE or WSS) around the
% optimal parameter values in the combinations qs1-qs2, b01-delU1, and
% b02-delU2.
%
% Last modified:
% - 2021-03-17, HA: Added comments
% - 2021-03-12, HA: Initial creation
%
% Input arguments:
% - x, y, z:             Pressure (x), Temperature (y), and adsorbed amount
%                        (z) data from experiments
%
% - nbins:               Number of groups data is binned to based on
%                        pressure (x) range in WSS method
%
% - isothermModel:       Isotherm model for which the data needs to be
%                        fitted ('DSL' or 'DSS')
%
% - fittingMethod:       Fitting method used for objective function
%                        generation (MLE or WSS)
%
% Output arguments:
% - Z:                   Matrices of objective function values for 3
%                        combinations of parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Z = generateObjfunContour(x,y,z,nbins,isothermModel,fittingMethod, isoRef,parVals)
switch isothermModel
    % For DSL model
    case 'DSL'
        % Determine number of points to consider in each axis
        nPoints = 200;
        % Obtain optimal parameter values from input
        qs1 = parVals(1)./isoRef(1);
        qs2 = parVals(2)./isoRef(2);
        b01 =  parVals(3)./isoRef(3);
        b02 =  parVals(4)./isoRef(4);
        delU1 =  parVals(5)./isoRef(5);
        delU2= parVals(6)./isoRef(6);
        % Define lower and upper bounds for parameters for contour plots
        parLB = 0.1.*parVals.*isoRef;
        parUB = 1.9.*parVals.*isoRef;
        % Create matrix for parameter range
        parRange = [];
        for jj = 1:length(parVals)
            parRange(:,jj) = linspace(parLB(jj),parUB(jj),nPoints);
        end
        % Obtain individual parameter ranges
        qs1vals = parRange(:,1);
        qs2vals = parRange(:,2);
        b01vals = parRange(:,3);
        b02vals = parRange(:,4);
        delU1vals = parRange(:,5);
        delU2vals = parRange(:,6);
        % create matrix for input and output mesh grids
        Z=zeros(nPoints,nPoints,1,3);
        X=zeros(nPoints,nPoints,1,3);
        Y=zeros(nPoints,nPoints,1,3);
        %% qs1 vs qs2
        % obtain objective function values for qs1 and qs2 parameter variation from
        % optima for WSS method
        [X(:,:,1,1),Y(:,:,1,1)] = meshgrid(qs1vals,qs2vals);
        for jj = 1:nPoints
            for kk =1:nPoints
                switch fittingMethod
                    case 'WSS'
                        Z(jj,kk,1,1) = generateWSSfun(x,y,z, nbins,  isothermModel, isoRef, qs1vals(jj), qs2vals(kk), b01, b02, delU1, delU2);
                    case 'MLE'
                        nbins = 1;
                        Z(jj,kk,1,1) = generateMLEfun(x,y,z, nbins,  isothermModel, isoRef, qs1vals(jj), qs2vals(kk), b01, b02, delU1, delU2);
                end
            end
        end
        
        %% b01 vs delU1
        % obtain objective function values for b01 and delU1 parameter variation from
        % optima for WSS method
        [X(:,:,1,2),Y(:,:,1,2)] = meshgrid(b01vals,delU1vals);
        for jj = 1:nPoints
            for kk =1:nPoints
                switch fittingMethod
                    case 'WSS'
                        Z(jj,kk,1,2) = generateWSSfun(x,y,z, nbins,  isothermModel, isoRef, qs1, qs2, b01vals(jj), b02, delU1vals(kk), delU2);
                    case 'MLE'
                        nbins = 1;
                        Z(jj,kk,1,2) = generateMLEfun(x,y,z, nbins,  isothermModel, isoRef, qs1, qs2, b01vals(jj), b02, delU1vals(kk), delU2);
                end
            end
        end
        
        %% b02 vs delU2
        % obtain objective function values for b02 and delU2 parameter variation from
        % optima for WSS method
        [X(:,:,1,3),Y(:,:,1,3)] = meshgrid(b02vals,delU2vals);
        for jj = 1:nPoints
            for kk =1:nPoints
                switch fittingMethod
                    case 'WSS'
                        Z(jj,kk,1,3) = generateWSSfun(x,y,z, nbins,  isothermModel, isoRef, qs1, qs2, b01, b02vals(jj), delU1, delU2vals(kk));
                    case 'MLE'
                        nbins = 1;
                        Z(jj,kk,1,3) = generateMLEfun(x,y,z, nbins,  isothermModel, isoRef, qs1, qs2, b01, b02vals(jj), delU1, delU2vals(kk));
                end
            end
        end
        
        % Plot figure with the contour subplots
        figure
        for mm = 1:3
            subplot(1,3,mm)
            contourf(X(:,:,1,mm),Y(:,:,1,mm),Z(:,:,1,mm),50,'--','LineWidth',0.2);
            box on
            set(gca,'YScale','linear','XScale','linear','FontSize',10,'LineWidth',1)
            grid on; axis square
            colorbar('southoutside');
            hold on
            switch mm
                case 1
                    plot(qs1,qs2,'r*');
                    xlabel('qs1 [mol/kg]');
                    ylabel('qs2 [mol/kg]');
                case 2
                    plot(b01,delU1,'r*');
                    xlabel('b01 [1/bar]');
                    ylabel('delU1 [J/mol]');
                    title(['Dual-site Langmuir with ',fittingMethod])
                case 3
                    plot(b02,delU2,'r*');
                    xlabel('b02 [1/bar]');
                    ylabel('delU2 [J/mol]');
            end
        end
        
        % For DSL model
    case 'DSS'
        nPoints = 200;
        qs1 = parVals(1)./isoRef(1);
        qs2 = parVals(2)./isoRef(2);
        b01 =  parVals(3)./isoRef(3);
        b02 =  parVals(4)./isoRef(4);
        delU1 =  parVals(5)./isoRef(5);
        delU2= parVals(6)./isoRef(6);
        gamma= parVals(6)./isoRef(7);
        % Define lower and upper bounds for parameters for contour plots
        parLB = 0.1.*parVals.*isoRef;
        parUB = 1.9.*parVals.*isoRef;
        
        parRange = [];
        for jj = 1:length(parVals)
            parRange(:,jj) = linspace(parLB(jj),parUB(jj),nPoints);
        end
        
        qs1vals = parRange(:,1);
        qs2vals = parRange(:,2);
        b01vals = parRange(:,3);
        b02vals = parRange(:,4);
        delU1vals = parRange(:,5);
        delU2vals = parRange(:,6);
        gammavals = parRange(:,7);
        
        Z=zeros(nPoints,nPoints,1,3);
        X=zeros(nPoints,nPoints,1,3);
        Y=zeros(nPoints,nPoints,1,3);
        
        %% qs1 vs qs2
        [X(:,:,1,1),Y(:,:,1,1)] = meshgrid(qs1vals,qs2vals);
        for jj = 1:nPoints
            for kk =1:nPoints
                switch fittingMethod
                    case 'WSS'
                        Z(jj,kk,1,1) = generateWSSfun(x,y,z, nbins,  isothermModel, isoRef, qs1vals(jj), qs2vals(kk), b01, b02, delU1, delU2, gamma);
                    case 'MLE'
                        nbins = 1;
                        Z(jj,kk,1,1) = generateMLEfun(x,y,z, nbins,  isothermModel, isoRef, qs1vals(jj), qs2vals(kk), b01, b02, delU1, delU2, gamma);
                end
            end
        end
        
        %% b01 vs delU1
        [X(:,:,1,2),Y(:,:,1,2)] = meshgrid(b01vals,delU1vals);
        for jj = 1:nPoints
            for kk =1:nPoints
                switch fittingMethod
                    case 'WSS'
                        Z(jj,kk,1,2) = generateWSSfun(x,y,z, nbins,  isothermModel, isoRef, qs1, qs2, b01vals(jj), b02, delU1vals(kk), delU2, gamma);
                    case 'MLE'
                        nbins = 1;
                        Z(jj,kk,1,2) = generateMLEfun(x,y,z, nbins,  isothermModel, isoRef, qs1, qs2, b01vals(jj), b02, delU1vals(kk), delU2, gamma);
                end
            end
        end
        
        %% b02 vs delU2
        [X(:,:,1,3),Y(:,:,1,3)] = meshgrid(b02vals,delU2vals);
        for jj = 1:nPoints
            for kk =1:nPoints
                switch fittingMethod
                    case 'WSS'
                        Z(jj,kk,1,3) = generateWSSfun(x,y,z, nbins,  isothermModel, isoRef, qs1, qs2, b01, b02vals(jj), delU1, delU2vals(kk), gamma);
                    case 'MLE'
                        nbins = 1;
                        Z(jj,kk,1,3) = generateMLEfun(x,y,z, nbins,  isothermModel, isoRef, qs1, qs2, b01, b02vals(jj), delU1, delU2vals(kk), gamma);
                end
            end
        end
        
        figure
        for mm = 1:3
            subplot(1,3,mm)
            contourf(X(:,:,1,mm),Y(:,:,1,mm),Z(:,:,1,mm),50,'--','LineWidth',0.2);
            colorbar('southoutside');
            hold on
            box on
            set(gca,'YScale','linear','XScale','linear','FontSize',10,'LineWidth',1)
            grid on; axis square
            switch mm
                case 1
                    plot(qs1,qs2,'r*');
                    xlabel('qs1 [mol/kg]');
                    ylabel('qs2 [mol/kg]');
                case 2
                    plot(b01,delU1,'r*');
                    xlabel('b01 [1/bar]');
                    ylabel('delU1 [J/mol]');
                    title(['Dual-site Sips with ',fittingMethod])
                case 3
                    plot(b02,delU2,'r*');
                    xlabel('b02 [1/bar]');
                    ylabel('delU2 [J/mol]');
            end
        end
end
end