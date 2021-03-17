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
% Obtain the uncertainty bounds for the isotherm using the isotherm
% parameter confidence intervals
%
% Last modified:
% - 2021-03-17, HA: Added comments
% - 2021-03-15, HA: Initial creation
%
% Input arguments:
% - x, y:                Pressure (x), and Temperature (y) from experiments
%
% - isothermModel:       Isotherm model for which the data needs to be
%                        fitted ('DSL' or 'DSS')
%
% - parVals:             Optimal fitting parameter values for the isotherms
%
% - conRange95:          Confidence intervals for the parameters
%
% Output arguments:
% - outScatter:          Matrix of points for uncertainty bound plots
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [outScatter]=generateUncertaintySpread(x,y,isothermModel,parVals,conRange95)
% Decide number of samplint points for q at each pressure point
nPoints = 30;
% Decide the range of pressure for the uncertainty spread calculation
Pvals=linspace(0,max(x),400);
% Obtain a vector of the temperatures present in the input data
Tvals=unique(y);
% Calculate the uncertainty spread for q in the pressure range
switch isothermModel
    % for the dual-site langmuir model
    case 'DSL'
        % obtain optimal parameter values from the input
        qs1 = parVals(1);
        qs2 = parVals(2);
        b01 =  parVals(3);
        b02 =  parVals(4);
        delU1 = parVals(5);
        delU2= parVals(6);
        % create 3 dimensional array for the output data
        qeqUnc=zeros(3,nPoints*length(Pvals),length(Tvals));
        % populate the first and second rows of the output array with the
        % temperature and pressure range
        for mm = 1:length(Tvals)
            for jj = 1:length(Pvals)
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj+1))
                    qeqUnc(1,kk,mm)=Pvals(jj);
                    qeqUnc(2,kk,mm)=Tvals(mm);
                end
            end
        end
        % Generate a random set of numbers between -1 and 1 using
        % latin-hypercube sampling in a matrix with nPoints rows and a
        % column for each variable
        lhsMat = 2*lhsdesign(nPoints,6)-1;
        % Populate the third row of the output array with the q values
        % calculated using the the values of parameters within their
        % respective uncertainty bounds
        for mm = 1:length(Tvals)
            for jj = 1:length(Pvals)
                hh = 0;
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj))
                    hh=hh+1;
                    % Determine values for each parameter within it's
                    % respective bounds
                    qs1_unc = qs1+conRange95(1)*lhsMat(hh,1);
                    qs2_unc = qs2+conRange95(2)*lhsMat(hh,2);
                    b01_unc = b01+conRange95(3)*lhsMat(hh,3);
                    b02_unc = b02+conRange95(4)*lhsMat(hh,4);
                    delU1_unc = delU1+conRange95(5)*lhsMat(hh,5);
                    delU2_unc = delU2+conRange95(6)*lhsMat(hh,6);
                    % Obtain P and T values corresponding to this data
                    % point
                    P = qeqUnc(1,kk,mm);
                    T = qeqUnc(2,kk,mm);
                    % Calculate q corresponding to the parameter values
                    % obtained above
                    qeqUnc(3,kk,mm)=qs1_unc.*(b01_unc.*P.*exp(delU1_unc./(8.314.*T)))./(1+(b01_unc.*P.*exp(delU1_unc./(8.314.*T)))) ...
                        + qs2_unc.*(b02_unc.*P.*exp(delU2_unc./(8.314.*T)))./(1+(b02_unc.*P.*exp(delU2_unc./(8.314.*T))));
                end
            end
        end
        % Obtain the output in a matrix [nPoints*length(Pvals) x 3]
        outScatter=[qeqUnc(1,:);qeqUnc(3,:); qeqUnc(2,:)];
        % Transpose of the output
        outScatter=outScatter';
        
    case 'DSS'
        qs1 = parVals(1);
        qs2 = parVals(2);
        b01 =  parVals(3);
        b02 =  parVals(4);
        delU1 =  parVals(5);
        delU2= parVals(6);
        gamma= parVals(7);
        
        qeqUnc=zeros(3,nPoints*length(Pvals),length(Tvals));
        for mm = 1:length(Tvals)
            for jj = 1:length(Pvals)
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj+1))
                    qeqUnc(1,kk,mm)=Pvals(jj);
                    qeqUnc(2,kk,mm)=Tvals(mm);
                end
            end
        end
        
        lhsMat = 2*lhsdesign(nPoints,7)-1;
        
        for mm = 1:length(Tvals)
            for jj = 1:length(Pvals)
                hh = 0;
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj))
                    hh=hh+1;
                    qs1_unc = qs1+conRange95(1)*lhsMat(hh,1);
                    qs2_unc = qs2+conRange95(2)*lhsMat(hh,2);
                    b01_unc = b01+conRange95(3)*lhsMat(hh,3);
                    b02_unc = b02+conRange95(4)*lhsMat(hh,4);
                    delU1_unc = delU1+conRange95(5)*lhsMat(hh,5);
                    delU2_unc = delU2+conRange95(6)*lhsMat(hh,6);
                    gamma_unc = gamma+conRange95(7)*lhsMat(hh,7);
                    P = qeqUnc(1,kk,mm);
                    T = qeqUnc(2,kk,mm);
                    qeqUnc(3,kk,mm)=qs1_unc.*(b01_unc.*P.*exp(delU1_unc./(8.314.*T))).^gamma_unc./(1+(b01_unc.*P.*exp(delU1_unc./(8.314.*T))).^gamma_unc) ...
                        + qs2_unc.*(b02_unc.*P.*exp(delU2_unc./(8.314.*T))).^gamma_unc./(1+(b02_unc.*P.*exp(delU2_unc./(8.314.*T))).^gamma_unc);
                end
            end
        end
        
        outScatter=[qeqUnc(1,:);qeqUnc(3,:); qeqUnc(2,:)];
        outScatter=outScatter';
end
end