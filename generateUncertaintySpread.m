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
% - 2021-10-04, HA: Add extremes of the uncertainty spread as an output
% - 2021-03-17, HA: Add comments
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
function [outScatter,uncBounds]=generateUncertaintySpread(x,y,z,isothermModel,parVals,conRange95,varargin)
% Decide number of samplint points for q at each pressure point
nPoints = 100;
% Decide the range of pressure for the uncertainty spread calculation
Pvals = linspace(0,max(x)*1.5,500);
% Decide the range of loading for the uncertainty spread calculation
% (virial only)
qvals = linspace(0,max(z)*1.5,8000);
% Calculate the uncertainty spread for q in the pressure range
switch isothermModel
    % for the dual-site langmuir model
    case 'STATZ'
        % Obtain a vector of the temperatures present in the input data
        Tvals = unique(y);
        uncBounds = [];
        uncBounds(:,1) = repmat(Pvals',length(Tvals),1);
        
        % obtain optimal parameter values from the input
        omega = parVals(1);
        beta = parVals(2);
        b01 = parVals(3);
        delU1 = parVals(4);
        vc = varargin{1};
        
        % create 3 dimensional array for the output data
        qeqUnc=zeros(3,nPoints*length(Pvals),length(Tvals));
        % populate the first and second rows of the output array with the
        % temperature and pressure range
        for mm = 1:length(Tvals)
            for jj = 1:length(Pvals)
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj+1))
                    qeqUnc(1,kk,mm) = Pvals(jj);
                    qeqUnc(2,kk,mm) = Tvals(mm);
                end
            end
        end
        % Generate a random set of numbers between -1 and 1 using
        % latin-hypercube sampling in a matrix with nPoints rows and a
        % column for each variable
        lhsMat = 2*lhsdesign(nPoints,4)-1;
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
                    omega_unc = omega+conRange95(1)*lhsMat(hh,1);
                    beta_unc = beta+conRange95(2)*lhsMat(hh,2);
                    b01_unc = b01+conRange95(3)*lhsMat(hh,3);
                    delU1_unc = delU1+conRange95(4)*lhsMat(hh,4);
                    % Obtain P and T values corresponding to this data
                    % point
                    P = qeqUnc(1,kk,mm);
                    T = qeqUnc(2,kk,mm);
                    % Calculate q corresponding to the parameter values
                    % obtained above
                    qeqUnc(3,kk,mm) = computeStatZLoading(P,T,b01_unc,delU1_unc,beta_unc,omega_unc,vc);
                end
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,2)= min(qeqUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,3)= max(qeqUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,4)= Tvals(mm);
            end
        end
        % Obtain the output in a matrix [nPoints*length(Pvals) x 3]
        outScatter=[qeqUnc(1,:);qeqUnc(3,:); qeqUnc(2,:)];
        % Transpose of the output
        outScatter=outScatter';
    case 'STATZGATE'
        % Obtain a vector of the temperatures present in the input data
        Tvals = unique(y);
        uncBounds = [];
        uncBounds(:,1) = repmat(Pvals',length(Tvals),1);
        
        % obtain optimal parameter values from the input
        omega = parVals(1);
        beta = parVals(2);
        b01 = parVals(3);
        delU1 = parVals(4);
        kgate = parVals(5);
        cgate = parVals(6);
        gamma = parVals(7);
        vc1 = varargin{1};
        vc2 = varargin{2};
      
        % create 3 dimensional array for the output data
        qeqUnc=zeros(3,nPoints*length(Pvals),length(Tvals));
        % populate the first and second rows of the output array with the
        % temperature and pressure range
        for mm = 1:length(Tvals)
            for jj = 1:length(Pvals)
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj+1))
                    qeqUnc(1,kk,mm) = Pvals(jj);
                    qeqUnc(2,kk,mm) = Tvals(mm);
                end
            end
        end
        % Generate a random set of numbers between -1 and 1 using
        % latin-hypercube sampling in a matrix with nPoints rows and a
        % column for each variable
        lhsMat = 2*lhsdesign(nPoints,7)-1;
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
                    omega_unc = omega;
                    beta_unc = beta+conRange95(2)*lhsMat(hh,2);
                    b01_unc = b01+conRange95(3)*lhsMat(hh,3);
                    delU1_unc = delU1+conRange95(4)*lhsMat(hh,4);
                    kgate_unc = kgate+conRange95(5)*lhsMat(hh,5);
                    cgate_unc = cgate+conRange95(6)*lhsMat(hh,6);
                    gamma_unc = gamma+conRange95(7)*lhsMat(hh,7);

                    % Obtain P and T values corresponding to this data
                    % point
                    P = qeqUnc(1,kk,mm);
                    T = qeqUnc(2,kk,mm);
                    % Calculate q corresponding to the parameter values
                    % obtained above
                    qeqUnc(3,kk,mm) = computeStatZGATELoading(P,T,b01_unc,delU1_unc,beta_unc,omega_unc,kgate_unc,cgate_unc,gamma_unc,vc1,vc2);
                end
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,2)= min(qeqUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,3)= max(qeqUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,4)= Tvals(mm);
            end
        end
        % Obtain the output in a matrix [nPoints*length(Pvals) x 3]
        outScatter=[qeqUnc(1,:);qeqUnc(3,:); qeqUnc(2,:)];
        % Transpose of the output
        outScatter=outScatter';
    case 'GAB'
        % Obtain a vector of the temperatures present in the input data
        Tvals = unique(y);
        uncBounds = [];
        uncBounds(:,1) = repmat(Pvals',length(Tvals),1);
        
        % obtain optimal parameter values from the input
        qs1 = parVals(1);
        parC = parVals(2);
        parD = parVals(3);
        parF = parVals(4);
        parG = parVals(5);
        
        % create 3 dimensional array for the output data
        qeqUnc=zeros(3,nPoints*length(Pvals),length(Tvals));
        % populate the first and second rows of the output array with the
        % temperature and pressure range
        for mm = 1:length(Tvals)
            for jj = 1:length(Pvals)
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj+1))
                    qeqUnc(1,kk,mm) = Pvals(jj);
                    qeqUnc(2,kk,mm) = Tvals(mm);
                end
            end
        end
        % Generate a random set of numbers between -1 and 1 using
        % latin-hypercube sampling in a matrix with nPoints rows and a
        % column for each variable
        lhsMat = 2*lhsdesign(nPoints,5)-1;
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
                    qs1_unc  = qs1+conRange95(1)*lhsMat(hh,1);
                    parC_unc = parC+conRange95(2)*lhsMat(hh,2);
                    parD_unc = parD+conRange95(3)*lhsMat(hh,3);
                    parF_unc = parF+conRange95(4)*lhsMat(hh,4);
                    parG_unc = parG+conRange95(5)*lhsMat(hh,5);
                    % Obtain P and T values corresponding to this data
                    % point
                    P = qeqUnc(1,kk,mm);
                    T = qeqUnc(2,kk,mm);
                    % Calculate q corresponding to the parameter values
                    % obtained above
                    qeqUnc(3,kk,mm) = computeGABLoading(P,T,qs1_unc,parC_unc,parD_unc,parF_unc,parG_unc);
                end
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,2)= min(qeqUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,3)= max(qeqUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,4)= Tvals(mm);
            end
        end
        % Obtain the output in a matrix [nPoints*length(Pvals) x 3]
        outScatter=[qeqUnc(1,:);qeqUnc(3,:); qeqUnc(2,:)];
        % Transpose of the output
        outScatter=outScatter';
    case 'UNIV6'
        % Obtain a vector of the temperatures present in the input data
        Tvals = unique(y);
        uncBounds = [];
        uncBounds(:,1) = repmat(Pvals',length(Tvals),1);
        
        % obtain optimal parameter values from the input
        qs = parVals(1);
        a1 = parVals(2);
        a2 = parVals(3);
        a3 = parVals(4);
        e01 = parVals(5);
        e02 = parVals(6);
        e03 = parVals(7);
        e04 = parVals(8);
        m1 = parVals(9);
        m2 = parVals(10);
        m3 = parVals(11);
        m4 = parVals(12);
        
        % create 3 dimensional array for the output data
        qeqUnc=zeros(3,nPoints*length(Pvals),length(Tvals));
        % populate the first and second rows of the output array with the
        % temperature and pressure range
        for mm = 1:length(Tvals)
            for jj = 1:length(Pvals)
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj+1))
                    qeqUnc(1,kk,mm) = Pvals(jj);
                    qeqUnc(2,kk,mm) = Tvals(mm);
                end
            end
        end
        % Generate a random set of numbers between -1 and 1 using
        % latin-hypercube sampling in a matrix with nPoints rows and a
        % column for each variable
        lhsMat = 2*lhsdesign(nPoints,12)-1;
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
                    qs_unc = qs+conRange95(1)*lhsMat(hh,1);
                    a1_unc = a1+conRange95(2)*lhsMat(hh,2);
                    a2_unc = a2+conRange95(3)*lhsMat(hh,3);
                    a3_unc = a3+conRange95(4)*lhsMat(hh,4);
                    e01_unc = e01+conRange95(5)*lhsMat(hh,5);
                    e02_unc = e02+conRange95(6)*lhsMat(hh,6);
                    e03_unc = e03+conRange95(7)*lhsMat(hh,7);
                    e04_unc = e04+conRange95(8)*lhsMat(hh,8);
                    m1_unc = m1+conRange95(9)*lhsMat(hh,9);
                    m2_unc = m2+conRange95(10)*lhsMat(hh,10);
                    m3_unc = m3+conRange95(11)*lhsMat(hh,11);
                    m4_unc = m4+conRange95(12)*lhsMat(hh,12);
                    % Obtain P and T values corresponding to this data
                    % point
                    P = qeqUnc(1,kk,mm);
                    T = qeqUnc(2,kk,mm);
                    % Calculate q corresponding to the parameter values
                    % obtained above
                    qeqUnc(3,kk,mm) = computeUNIV6Loading(P,T, qs_unc, a1_unc, a2_unc, a3_unc, e01_unc, e02_unc, e03_unc, e04_unc, m1_unc, m2_unc, m3_unc, m4_unc);
                end
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,2)= min(qeqUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,3)= max(qeqUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,4)= Tvals(mm);
            end
        end
        % Obtain the output in a matrix [nPoints*length(Pvals) x 3]
        outScatter=[qeqUnc(1,:);qeqUnc(3,:); qeqUnc(2,:)];
        % Transpose of the output
        outScatter=outScatter';
    case 'DSL'
        % Obtain a vector of the temperatures present in the input data
        Tvals = unique(y);
        uncBounds = [];
        uncBounds(:,1) = repmat(Pvals',length(Tvals),1);
        
        % obtain optimal parameter values from the input
        qs1 = parVals(1);
        qs2 = parVals(2);
        b01 =  parVals(3);
        b02 =  parVals(4);
        delU1 = parVals(5);
        delU2 = parVals(6);
        % create 3 dimensional array for the output data
        qeqUnc=zeros(3,nPoints*length(Pvals),length(Tvals));
        % populate the first and second rows of the output array with the
        % temperature and pressure range
        for mm = 1:length(Tvals)
            for jj = 1:length(Pvals)
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj+1))
                    qeqUnc(1,kk,mm) = Pvals(jj);
                    qeqUnc(2,kk,mm) = Tvals(mm);
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
                    qeqUnc(3,kk,mm) = qs1_unc.*(b01_unc.*P.*exp(delU1_unc./(8.314.*T)))./(1+(b01_unc.*P.*exp(delU1_unc./(8.314.*T)))) ...
                        + qs2_unc.*(b02_unc.*P.*exp(delU2_unc./(8.314.*T)))./(1+(b02_unc.*P.*exp(delU2_unc./(8.314.*T))));
                end
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,2)= min(qeqUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,3)= max(qeqUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,4)= Tvals(mm);
            end
        end
        % Obtain the output in a matrix [nPoints*length(Pvals) x 3]
        outScatter=[qeqUnc(1,:);qeqUnc(3,:); qeqUnc(2,:)];
        % Transpose of the output
        outScatter=outScatter';
    case 'DSL2'
        % Obtain a vector of the temperatures present in the input data
        Tvals = unique(y);
        uncBounds = [];
        uncBounds(:,1) = repmat(Pvals',length(Tvals),1);
        
        % obtain optimal parameter values from the input
        qs1a = parVals(1);
        qs2a = parVals(2);
        qs1b = parVals(3);
        qs2b = parVals(4);
        b01 =  parVals(5);
        b02 =  parVals(6);
        delU1 = parVals(7);
        delU2 = parVals(8);
        % create 3 dimensional array for the output data
        qeqUnc=zeros(3,nPoints*length(Pvals),length(Tvals));
        % populate the first and second rows of the output array with the
        % temperature and pressure range
        for mm = 1:length(Tvals)
            for jj = 1:length(Pvals)
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj+1))
                    qeqUnc(1,kk,mm) = Pvals(jj);
                    qeqUnc(2,kk,mm) = Tvals(mm);
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
            for jj = 1:length(Pvals)
                hh = 0;
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj))
                    hh=hh+1;
                    % Determine values for each parameter within it's
                    % respective bounds
                    qs1a_unc =   qs1a+conRange95(1)*lhsMat(hh,1);
                    qs2a_unc =   qs2a+conRange95(2)*lhsMat(hh,2);
                    qs1b_unc =   qs1b+conRange95(3)*lhsMat(hh,3);
                    qs2b_unc =   qs2b+conRange95(4)*lhsMat(hh,4);
                    b01_unc =     b01+conRange95(5)*lhsMat(hh,5);
                    b02_unc =     b02+conRange95(6)*lhsMat(hh,6);
                    delU1_unc = delU1+conRange95(7)*lhsMat(hh,7);
                    delU2_unc = delU2+conRange95(8)*lhsMat(hh,8);
                    % Obtain P and T values corresponding to this data
                    % point
                    P = qeqUnc(1,kk,mm);
                    T = qeqUnc(2,kk,mm);
                    % Calculate q corresponding to the parameter values
                    % obtained above
                    qs1_unc = qs1a_unc + qs1b_unc./T;
                    qs2_unc = qs2a_unc + qs2b_unc./T;
                    qeqUnc(3,kk,mm) = qs1_unc.*(b01_unc.*P.*exp(delU1_unc./(8.314.*T)))./(1+(b01_unc.*P.*exp(delU1_unc./(8.314.*T)))) ...
                        + qs2_unc.*(b02_unc.*P.*exp(delU2_unc./(8.314.*T)))./(1+(b02_unc.*P.*exp(delU2_unc./(8.314.*T))));
                end
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,2)= min(qeqUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,3)= max(qeqUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,4)= Tvals(mm);
            end
        end
        % Obtain the output in a matrix [nPoints*length(Pvals) x 3]
        outScatter=[qeqUnc(1,:);qeqUnc(3,:); qeqUnc(2,:)];
        % Transpose of the output
        outScatter=outScatter';
    case 'HDSL'
        % Obtain a vector of the temperatures present in the input data
        Tvals = unique(y);
        uncBounds = [];
        uncBounds(:,1) = repmat(Pvals',length(Tvals),1);
        
        % obtain optimal parameter values from the input
        qs1 = parVals(1);
        qs2 = parVals(2);
        b01 =  parVals(3);
        b02 =  parVals(4);
        delU1 = parVals(5);
        delU2 = parVals(6);
        qsH =  parVals(7);
        b0H = parVals(8);
        delUH = parVals(9);
        % create 3 dimensional array for the output data
        qeqUnc=zeros(3,nPoints*length(Pvals),length(Tvals));
        % populate the first and second rows of the output array with the
        % temperature and pressure range
        for mm = 1:length(Tvals)
            for jj = 1:length(Pvals)
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj+1))
                    qeqUnc(1,kk,mm) = Pvals(jj);
                    qeqUnc(2,kk,mm) = Tvals(mm);
                end
            end
        end
        % Generate a random set of numbers between -1 and 1 using
        % latin-hypercube sampling in a matrix with nPoints rows and a
        % column for each variable
        lhsMat = 2*lhsdesign(nPoints,9)-1;
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
                    qsH_unc = qsH+conRange95(7)*lhsMat(hh,7);
                    b0H_unc = b0H+conRange95(8)*lhsMat(hh,8);
                    delUH_unc = delUH+conRange95(9)*lhsMat(hh,9);
                    % Obtain P and T values corresponding to this data
                    % point
                    P = qeqUnc(1,kk,mm);
                    T = qeqUnc(2,kk,mm);
                    % Calculate q corresponding to the parameter values
                    % obtained above
                    qeqUnc(3,kk,mm) = qs1_unc.*(b01_unc.*P.*exp(delU1_unc./(8.314.*T)))./(1+(b01_unc.*P.*exp(delU1_unc./(8.314.*T)))) ...
                        + qs2_unc.*(b02_unc.*P.*exp(delU2_unc./(8.314.*T)))./(1+(b02_unc.*P.*exp(delU2_unc./(8.314.*T)))) ...
                        + qsH_unc.*(b0H_unc.*P.*exp(delUH_unc./(8.314.*T)));
                end
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,2)= min(qeqUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,3)= max(qeqUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,4)= Tvals(mm);
            end
        end
        % Obtain the output in a matrix [nPoints*length(Pvals) x 3]
        outScatter=[qeqUnc(1,:);qeqUnc(3,:); qeqUnc(2,:)];
        % Transpose of the output
        outScatter=outScatter';
    case 'TSL'
        % Obtain a vector of the temperatures present in the input data
        Tvals = unique(y);
        uncBounds = [];
        uncBounds(:,1) = repmat(Pvals',length(Tvals),1);
        
        % obtain optimal parameter values from the input
        qs1 = parVals(1);
        qs2 = parVals(2);
        b01 =  parVals(3);
        b02 =  parVals(4);
        delU1 = parVals(5);
        delU2 = parVals(6);
        qs3 =  parVals(7);
        b03 = parVals(8);
        delU3 = parVals(9);
        % create 3 dimensional array for the output data
        qeqUnc=zeros(3,nPoints*length(Pvals),length(Tvals));
        % populate the first and second rows of the output array with the
        % temperature and pressure range
        for mm = 1:length(Tvals)
            for jj = 1:length(Pvals)
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj+1))
                    qeqUnc(1,kk,mm) = Pvals(jj);
                    qeqUnc(2,kk,mm) = Tvals(mm);
                end
            end
        end
        % Generate a random set of numbers between -1 and 1 using
        % latin-hypercube sampling in a matrix with nPoints rows and a
        % column for each variable
        lhsMat = 2*lhsdesign(nPoints,9)-1;
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
                    qs3_unc = qs3+conRange95(7)*lhsMat(hh,7);
                    b03_unc = b03+conRange95(8)*lhsMat(hh,8);
                    delU3_unc = delU3+conRange95(9)*lhsMat(hh,9);
                    % Obtain P and T values corresponding to this data
                    % point
                    P = qeqUnc(1,kk,mm);
                    T = qeqUnc(2,kk,mm);
                    % Calculate q corresponding to the parameter values
                    % obtained above
                    qeqUnc(3,kk,mm) = qs1_unc.*(b01_unc.*P.*exp(delU1_unc./(8.314.*T)))./(1+(b01_unc.*P.*exp(delU1_unc./(8.314.*T)))) ...
                        + qs2_unc.*(b02_unc.*P.*exp(delU2_unc./(8.314.*T)))./(1+(b02_unc.*P.*exp(delU2_unc./(8.314.*T)))) ...
                        + qs3_unc.*(b03_unc.*P.*exp(delU3_unc./(8.314.*T)))./(1+(b03_unc.*P.*exp(delU3_unc./(8.314.*T))));
                end
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,2)= min(qeqUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,3)= max(qeqUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,4)= Tvals(mm);
            end
        end
        % Obtain the output in a matrix [nPoints*length(Pvals) x 3]
        outScatter=[qeqUnc(1,:);qeqUnc(3,:); qeqUnc(2,:)];
        % Transpose of the output
        outScatter=outScatter';
        % For sips model
    case 'DSS'
        % Obtain a vector of the temperatures present in the input data
        Tvals = unique(y);
        uncBounds = [];
        uncBounds(:,1) = repmat(Pvals',length(Tvals),1);
        
        qs1 = parVals(1);
        qs2 = parVals(2);
        b01 =  parVals(3);
        b02 =  parVals(4);
        delU1 = parVals(5);
        delU2 = parVals(6);
        gamma = parVals(7);
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
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,2)= min(qeqUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,3)= max(qeqUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,4)= Tvals(mm);
            end
        end
        outScatter=[qeqUnc(1,:);qeqUnc(3,:); qeqUnc(2,:)];
        outScatter=outScatter';
    case 'TOTH'
        % Obtain a vector of the temperatures present in the input data
        Tvals = unique(y);
        uncBounds = [];
        uncBounds(:,1) = repmat(Pvals',length(Tvals),1);
        
        qs1 = parVals(1);
        qs2 = parVals(2);
        b01 =  parVals(3);
        b02 =  parVals(4);
        delU1 = parVals(5);
        delU2 = parVals(6);
        toth = parVals(7);
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
                    toth_unc = toth+conRange95(7)*lhsMat(hh,7);
                    P = qeqUnc(1,kk,mm);
                    T = qeqUnc(2,kk,mm);
                    qeqUnc(3,kk,mm)=qs1_unc.*b01_unc.*P.*exp(delU1_unc./(8.314.*T))./(1+(b01_unc.*P.*exp(delU1_unc./(8.314.*T))).^toth_unc).^(1./toth_unc);
                end
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,2)= min(qeqUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,3)= max(qeqUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,4)= Tvals(mm);
            end
        end
        outScatter=[qeqUnc(1,:);qeqUnc(3,:); qeqUnc(2,:)];
        outScatter=outScatter';
    case 'TOTH2'
        % Obtain a vector of the temperatures present in the input data
        Tvals = unique(y);
        uncBounds = [];
        uncBounds(:,1) = repmat(Pvals',length(Tvals),1);
        
        qs1 = parVals(1);
        qs2 = parVals(2);
        b01 =  parVals(3);
        b02 =  parVals(4);
        delU1 = parVals(5);
        delU2 = parVals(6);
        toth0 = parVals(7);
        totha = parVals(8);
        qeqUnc=zeros(3,nPoints*length(Pvals),length(Tvals));
        for mm = 1:length(Tvals)
            for jj = 1:length(Pvals)
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj+1))
                    qeqUnc(1,kk,mm)=Pvals(jj);
                    qeqUnc(2,kk,mm)=Tvals(mm);
                end
            end
        end
        lhsMat = 2*lhsdesign(nPoints,8)-1;
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
                    toth0_unc = toth0+conRange95(7)*lhsMat(hh,7);
                    totha_unc = totha+conRange95(8)*lhsMat(hh,8);
                    P = qeqUnc(1,kk,mm);
                    T = qeqUnc(2,kk,mm);
                    qeqUnc(3,kk,mm)=qs1_unc.*b01_unc.*P.*exp(delU1_unc./(8.314.*T))./(1+(b01_unc.*P.*exp(delU1_unc./(8.314.*T))).^(toth0_unc + totha_unc.*(1-298./T))).^(1./(toth0_unc + totha_unc.*(1-298/T)));
                end
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,2)= min(qeqUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,3)= max(qeqUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,4)= Tvals(mm);
            end
        end
        outScatter=[qeqUnc(1,:);qeqUnc(3,:); qeqUnc(2,:)];
        outScatter=outScatter';
    case 'TOTH3'
        % Obtain a vector of the temperatures present in the input data
        Tvals = unique(y);
        uncBounds = [];
        uncBounds(:,1) = repmat(Pvals',length(Tvals),1);
        
        qs10 = parVals(1);
        qs2 = parVals(2);
        b01 =  parVals(3);
        b02 =  parVals(4);
        delU1 = parVals(5);
        delU2 = parVals(6);
        toth0 = parVals(7);
        totha = parVals(8);
        chi = parVals(9);
        qeqUnc=zeros(3,nPoints*length(Pvals),length(Tvals));
        for mm = 1:length(Tvals)
            for jj = 1:length(Pvals)
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj+1))
                    qeqUnc(1,kk,mm)=Pvals(jj);
                    qeqUnc(2,kk,mm)=Tvals(mm);
                end
            end
        end
        lhsMat = 2*lhsdesign(nPoints,9)-1;
        for mm = 1:length(Tvals)
            for jj = 1:length(Pvals)
                hh = 0;
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj))
                    hh=hh+1;
                    qs10_unc = qs10+conRange95(1)*lhsMat(hh,1);
                    qs2_unc = qs2+conRange95(2)*lhsMat(hh,2);
                    b01_unc = b01+conRange95(3)*lhsMat(hh,3);
                    b02_unc = b02+conRange95(4)*lhsMat(hh,4);
                    delU1_unc = delU1+conRange95(5)*lhsMat(hh,5);
                    delU2_unc = delU2+conRange95(6)*lhsMat(hh,6);
                    toth0_unc = toth0+conRange95(7)*lhsMat(hh,7);
                    totha_unc = totha+conRange95(8)*lhsMat(hh,8);
                    chi_unc = chi+conRange95(9)*lhsMat(hh,9);
                    P = qeqUnc(1,kk,mm);
                    T = qeqUnc(2,kk,mm);
                    qeqUnc(3,kk,mm)=qs10_unc.*exp(chi_unc*(1-T./298)).*b01_unc.*P.*exp(delU1_unc./(8.314.*T))./(1+(b01_unc.*P.*exp(delU1_unc./(8.314.*T))).^(toth0_unc + totha_unc.*(1-298./T))).^(1./(toth0_unc + totha_unc.*(1-298/T)));
                end
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,2)= min(qeqUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,3)= max(qeqUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,4)= Tvals(mm);
            end
        end
        outScatter=[qeqUnc(1,:);qeqUnc(3,:); qeqUnc(2,:)];
        outScatter=outScatter';
    case 'TOTHCHEM'
        % Obtain a vector of the temperatures present in the input data
        Tvals = unique(y);
        uncBounds = [];
        uncBounds(:,1) = repmat(Pvals',length(Tvals),1);
        
        qs10 = parVals(1);
        qs2 = parVals(2);
        b01 =  parVals(3);
        b02 =  parVals(4);
        delU1 = parVals(5);
        delU2 = parVals(6);
        toth0 = parVals(7);
        totha = parVals(8);
        chi = parVals(9);
        qsC = parVals(10);
        b0C = parVals(11);
        delUC = parVals(12);
        EaC = parVals(13);
        qeqUnc=zeros(3,nPoints*length(Pvals),length(Tvals));
        for mm = 1:length(Tvals)
            for jj = 1:length(Pvals)
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj+1))
                    qeqUnc(1,kk,mm)=Pvals(jj);
                    qeqUnc(2,kk,mm)=Tvals(mm);
                end
            end
        end
        lhsMat = 2*lhsdesign(nPoints,9)-1;
        for mm = 1:length(Tvals)
            for jj = 1:length(Pvals)
                hh = 0;
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj))
                    hh=hh+1;
                    qs10_unc = qs10+conRange95(1)*lhsMat(hh,1);
                    qs2_unc = qs2+conRange95(2)*lhsMat(hh,2);
                    b01_unc = b01+conRange95(3)*lhsMat(hh,3);
                    b02_unc = b02+conRange95(4)*lhsMat(hh,4);
                    delU1_unc = delU1+conRange95(5)*lhsMat(hh,5);
                    delU2_unc = delU2+conRange95(6)*lhsMat(hh,6);
                    toth0_unc = toth0+conRange95(7)*lhsMat(hh,7);
                    totha_unc = totha+conRange95(8)*lhsMat(hh,8);
                    chi_unc = chi+conRange95(9)*lhsMat(hh,9);
                    P = qeqUnc(1,kk,mm);
                    T = qeqUnc(2,kk,mm);
                    qeqUnc(3,kk,mm)=qs10_unc.*exp(chi_unc*(1-T./298)).*b01_unc.*P.*exp(delU1_unc./(8.314.*T))./(1+(b01_unc.*P.*exp(delU1_unc./(8.314.*T))).^(toth0_unc + totha_unc.*(1-298./T))).^(1./(toth0_unc + totha_unc.*(1-298/T))) ...
                        + exp(-EaC/(8.314.*T))*qsC.*b0C.*P.*exp(delUC./(8.314.*T))./(1+b0C.*P.*exp(delUC./(8.314.*T)));
                end
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,2)= min(qeqUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,3)= max(qeqUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(Pvals)*(mm-1)+(jj-1)+1,4)= Tvals(mm);
            end
        end
        outScatter=[qeqUnc(1,:);qeqUnc(3,:); qeqUnc(2,:)];
        outScatter=outScatter';
        % for the virial model
    case 'VIRIAL'
        % Obtain a vector of the temperatures present in the input data
        Tvals = unique(y);
        uncBounds = [];
        uncBounds(:,1) = repmat(qvals',length(Tvals),1);
        
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
        lnPUnc=zeros(3,nPoints*length(qvals),length(Tvals));
        % populate the first and second rows of the output array with the
        % temperature and pressure range
        for mm = 1:length(Tvals)
            for jj = 1:length(qvals)
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj+1))
                    lnPUnc(1,kk,mm) = qvals(jj);
                    lnPUnc(2,kk,mm) = Tvals(mm);
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
                    q = lnPUnc(1,kk,mm);
                    T = lnPUnc(2,kk,mm);
                    % Calculate q corresponding to the parameter values
                    % obtained above
                    lnPUnc(3,kk,mm) = log(q) + 1./T.*(a0_unc + a1_unc.*q + ...
                        a2_unc.*q^2 + a3_unc.*q^3) + b0_unc ...
                        + b1_unc.*q + b2_unc.*q.^2 + b3_unc.*q.^3;
                end
                uncBounds(length(qvals)*(mm-1)+(jj-1)+1,2)= min(lnPUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(qvals)*(mm-1)+(jj-1)+1,3)= max(lnPUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(qvals)*(mm-1)+(jj-1)+1,4)= Tvals(mm);
            end
        end
        % Obtain the output in a matrix [nPoints*length(Pvals) x 3]
        outScatter=[lnPUnc(1,:);lnPUnc(3,:); lnPUnc(2,:)];
        % Transpose of the output
        outScatter=outScatter';
    case 'VIRIAL2'
        % Obtain a vector of the temperatures present in the input data
        Tvals = unique(y);
        uncBounds = [];
        uncBounds(:,1) = repmat(qvals',length(Tvals),1);
        
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
        lnPUnc=zeros(3,nPoints*length(qvals),length(Tvals));
        % populate the first and second rows of the output array with the
        % temperature and pressure range
        for mm = 1:length(Tvals)
            for jj = 1:length(qvals)
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj+1))
                    lnPUnc(1,kk,mm) = qvals(jj);
                    lnPUnc(2,kk,mm) = Tvals(mm);
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
                    q = lnPUnc(1,kk,mm);
                    T = lnPUnc(2,kk,mm);
                    % Calculate q corresponding to the parameter values
                    % obtained above
                    lnPUnc(3,kk,mm) = log(q) + 1./T.*(a0_unc + a1_unc.*q + ...
                        a2_unc.*q^2 + a3_unc.*q^3) + b0_unc ...
                        + b1_unc.*q + b2_unc.*q.^2 + b3_unc.*q.^3;
                end
                uncBounds(length(qvals)*(mm-1)+(jj-1)+1,2)= min(lnPUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(qvals)*(mm-1)+(jj-1)+1,3)= max(lnPUnc(3,(nPoints*(jj-1)+1):(nPoints*(jj)),mm));
                uncBounds(length(qvals)*(mm-1)+(jj-1)+1,4)= Tvals(mm);
            end
        end
        % Obtain the output in a matrix [nPoints*length(Pvals) x 3]
        outScatter=[lnPUnc(1,:);lnPUnc(3,:); lnPUnc(2,:)];
        % Transpose of the output
        outScatter=outScatter';
end
end