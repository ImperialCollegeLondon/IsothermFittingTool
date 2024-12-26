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
% Obtain the uncertainty bounds for the heat of adsorption obtained using
% virial isotherm fitting
%
% Last modified:
% - 2023-03-23, HA: Initial creation
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
function [delHoutScatter,delHuncBounds]=generateUncertaintydelH(x,y,z,isothermModel,parVals,conRange95,varargin)
% Decide number of samplint points for q at each pressure point
nPoints = 500;
% Decide the range of pressure for the uncertainty spread calculation
Pvals = linspace(0,max(x),100);
% Decide the range of loading for the uncertainty spread calculation
% (virial only)
qvals = linspace(0,max(z),100);

switch isothermModel
 % for the virial model
    case 'VIRIAL'
        % Obtain a vector of the temperatures present in the input data
        Tvals = unique(y);
        Tvals = Tvals(1);
        delHuncBounds = [];
        delHuncBounds(:,1) = repmat(qvals',length(Tvals),1);
        
        % obtain optimal parameter values from the input
        a0 = parVals(1);
        a1 = parVals(2);
        a2 = parVals(3);
        a3 = parVals(4);
        b0 = parVals(5);
        b1 = parVals(6);
        b2 =0;
        b3 =0;
        % create 3 dimensional array for the output data
        delHUnc=zeros(3,nPoints*length(qvals),length(Tvals));
        % populate the first and second rows of the output array with the
        % temperature and pressure range
        for mm = 1:length(Tvals)
            for jj = 1:length(qvals)
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj+1))
                    delHUnc(1,kk,mm) = qvals(jj);
                    delHUnc(2,kk,mm) = Tvals(mm);
                end
            end
        end
        % Generate a random set of numbers between -1 and 1 using
        % latin-hypercube sampling in a matrix with nPoints rows and a
        % column for each variable
        lhsMat = 2*lhsdesign(nPoints,8)-1;
        % Populate the third row of the output array with the q values
        % calculated using the the values of parameters within their
        % respective uncertainty bounds
        for mm = 1:length(Tvals)
            for jj = 1:length(qvals)
                hh = 0;
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj))
                    hh=hh+1;
                    % Determine values for each parameter within it's
                    % respective bounds
                    a0_unc = a0+conRange95(1)*lhsMat(hh,1);
                    a1_unc = a1+conRange95(2)*lhsMat(hh,2);
                    a2_unc = a2+conRange95(3)*lhsMat(hh,3);
                    a3_unc = a3+conRange95(4)*lhsMat(hh,4);
                    b0_unc = b0+conRange95(5)*lhsMat(hh,5);
                    b1_unc = b1+conRange95(6)*lhsMat(hh,6);
                    b2_unc = b2+conRange95(7)*lhsMat(hh,7);
                    b3_unc = b3+conRange95(8)*lhsMat(hh,8);
                    % Obtain q and T values corresponding to this data
                    % point
                    q = delHUnc(1,kk,mm);
                    T = delHUnc(2,kk,mm);
                    % Calculate q corresponding to the parameter values
                    % obtained above
                    delHUnc(3,kk,mm) = -8.314.*(a0_unc + a1_unc.*q + ...
                        a2_unc.*q^2 + a3_unc.*q^3);
                end
                delHuncBounds(length(qvals)*(mm-1)+(jj-1)+1,2)= min(delHUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                delHuncBounds(length(qvals)*(mm-1)+(jj-1)+1,3)= max(delHUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                delHuncBounds(length(qvals)*(mm-1)+(jj-1)+1,4)= Tvals(mm);
            end
        end
        % Obtain the output in a matrix [nPoints*length(Pvals) x 3]
        delHoutScatter=[delHUnc(1,:);delHUnc(3,:); delHUnc(2,:)];
        % Transpose of the output
        delHoutScatter=delHoutScatter';
    case 'VIRIAL2'
        % Obtain a vector of the temperatures present in the input data
        Tvals = unique(y);
        Tvals = Tvals(1);
        delHuncBounds = [];
        delHuncBounds(:,1) = repmat(qvals',length(Tvals),1);
        
        % obtain optimal parameter values from the input
        a0 = parVals(1);
        a1 = parVals(2);
        a2 = parVals(3);
        a3 = parVals(4);
        b0 = parVals(5);
        b1 = parVals(6);
        b2 = parVals(7);
        b3 = parVals(8);
        % create 3 dimensional array for the output data
        delHUnc=zeros(3,nPoints*length(qvals),length(Tvals));
        % populate the first and second rows of the output array with the
        % temperature and pressure range
        for mm = 1:length(Tvals)
            for jj = 1:length(qvals)
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj+1))
                    delHUnc(1,kk,mm) = qvals(jj);
                    delHUnc(2,kk,mm) = Tvals(mm);
                end
            end
        end
        % Generate a random set of numbers between -1 and 1 using
        % latin-hypercube sampling in a matrix with nPoints rows and a
        % column for each variable
        lhsMat = 2*lhsdesign(nPoints,8)-1;
        % Populate the third row of the output array with the q values
        % calculated using the the values of parameters within their
        % respective uncertainty bounds
        for mm = 1:length(Tvals)
            for jj = 1:length(qvals)
                hh = 0;
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj))
                    hh=hh+1;
                    % Determine values for each parameter within it's
                    % respective bounds
                    a0_unc = a0+conRange95(1)*lhsMat(hh,1);
                    a1_unc = a1+conRange95(2)*lhsMat(hh,2);
                    a2_unc = a2+conRange95(3)*lhsMat(hh,3);
                    a3_unc = a3+conRange95(4)*lhsMat(hh,4);
                    b0_unc = b0+conRange95(5)*lhsMat(hh,5);
                    b1_unc = b1+conRange95(6)*lhsMat(hh,6);
                    b2_unc = b2+conRange95(7)*lhsMat(hh,7);
                    b3_unc = b3+conRange95(8)*lhsMat(hh,8);
                    % Obtain q and T values corresponding to this data
                    % point
                    q = delHUnc(1,kk,mm);
                    T = delHUnc(2,kk,mm);
                    % Calculate q corresponding to the parameter values
                    % obtained above
                    delHUnc(3,kk,mm) = -8.314.*(a0_unc + a1_unc.*q + ...
                        a2_unc.*q^2 + a3_unc.*q^3);
                end
                delHuncBounds(length(qvals)*(mm-1)+(jj-1)+1,2)= min(delHUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                delHuncBounds(length(qvals)*(mm-1)+(jj-1)+1,3)= max(delHUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                delHuncBounds(length(qvals)*(mm-1)+(jj-1)+1,4)= Tvals(mm);
            end
        end
        % Obtain the output in a matrix [nPoints*length(Pvals) x 3]
        delHoutScatter=[delHUnc(1,:);delHUnc(3,:); delHUnc(2,:)];
        % Transpose of the output
        delHoutScatter=delHoutScatter';

end