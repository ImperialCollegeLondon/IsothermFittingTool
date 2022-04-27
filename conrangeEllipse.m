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

function [conRange95] = conrangeEllipse(x,y,z,fitVals,fittingMethod,isoRef,isothermModel,varargin)
switch isothermModel
    % For DSL model
    case 'DSL'
        % Number of parameters
        Np = length(cell2mat(varargin(cell2mat(varargin)~=0)));
        % Calculate standard deviation of the data (not needed)
        stDevData = sqrt(1/(length(x)-length(Np)) * sum((z-fitVals).^2));
        % degree of variation of parameter to calculate sensitivity
        del = 0.000001;
        % Generate and solve global optimisation problem for confidence regions
        % based on isotherm model
        qs1 = varargin{1};
        qs2 = varargin{2};
        b01 = varargin{3};
        b02 = varargin{4};
        delU1 = varargin{5};
        delU2 = varargin{6};
        switch fittingMethod
            case 'WSS'
                % Create empty sensitivity matrix
                sensitivityMatrix = zeros(length(x),6);
                % Calculate sensitivity at every data point for each parameter
                for jj = 1:6
                    for kk = 1:length(x)
                        if jj == 1
                            sensitivityMatrix(kk,jj) = ((((1+del)*qs1)*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + qs2*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))))) - fitVals(kk))/(del*qs1);
                        elseif jj == 2
                            sensitivityMatrix(kk,jj) = ((qs1*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + ((1+del)*qs2)*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))))) - fitVals(kk))/(del*qs2);
                        elseif jj == 3
                            sensitivityMatrix(kk,jj) = ((qs1*(((1+del)*b01)*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(((1+del)*b01)*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + qs2*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))))) - fitVals(kk))/(del*b01);
                        elseif jj == 4
                            sensitivityMatrix(kk,jj) = ((qs1*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + qs2*(((1+del)*b02)*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(((1+del)*b02)*x(kk)*exp(delU2/(8.314*y(kk)))))) - fitVals(kk))/(del*b02);
                        elseif jj == 5
                            sensitivityMatrix(kk,jj) = ((qs1*(b01*x(kk)*exp((1+del)*delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp((1+del)*delU1/(8.314*y(kk))))) ...
                                + qs2*((b02)*x(kk)*exp(delU2/(8.314*y(kk))))/(1+((b02)*x(kk)*exp(delU2/(8.314*y(kk)))))) - fitVals(kk))/(del*delU1);
                        else
                            sensitivityMatrix(kk,jj) = ((qs1*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + qs2*((b02)*x(kk)*exp((1+del)*delU2/(8.314*y(kk))))/(1+((b02)*x(kk)*exp((1+del)*delU2/(8.314*y(kk)))))) - fitVals(kk))/(del*delU2);
                        end
                    end
                end
                % estimated Hessian Matrix for the data set (Non-linear parameter estimation
                % by Yonathan Bard (1974) pg. 178
                hessianMatrix = 1/stDevData^2*transpose(sensitivityMatrix)*sensitivityMatrix;
                % Confidence range given by chi squared distribution at Np degrees
                % of freedom (independent parameter conf intervals)
                conRange95 = sqrt(chi2inv(0.95,Np)./diag(hessianMatrix));
            case 'MLE'
                Np = 6;
                Nt = length(x);
                del = 1e-6;
                dlogMLE = [];
                d2logMLE = [];
                deltaplus1mat = eye(Np).*(del);
                deltamat = eye(Np).*(del);
                partemp = [qs1./isoRef(1), qs2./isoRef(2),b01./isoRef(3),b02./isoRef(4),delU1./isoRef(5),delU2./isoRef(6)];
                logMLE = @(par) -generateMLEfun(x, y, z, 1, 'DSL', isoRef, par(1), par(2), par(3), ...
                    par(4), par(5), par(6));
                
                for jj = 1:Np
                    for kk = 1:Np
                        partempnumj = partemp.*(1+deltaplus1mat(jj,:));
                        partempdenj = partemp.*deltamat(jj,:);
                        partempnumk = partemp.*(1+deltaplus1mat(kk,:));
                        partempdenk = partemp.*deltamat(kk,:);
                        partempnumjk = partemp.*(1+deltaplus1mat(jj,:) + deltaplus1mat(kk,:));
                        % Compute second derivative of logL for jj and kk
                        d2logMLE(jj,kk) = ((logMLE(partempnumjk)-logMLE(partempnumk))-(logMLE(partempnumj)-logMLE(partemp)))./(partempdenj(jj).*isoRef(jj).*partempdenk(kk).*isoRef(kk));
                    end
                end
                % estimated Hessian Matrix for the data set (Non-linear parameter estimation
                % by Yonathan Bard (1974) pg. 178 (Eqn 7-5-17)
                hessianMatrix =  -d2logMLE;
                conRange95 = sqrt(chi2inv(0.95,length(cell2mat(varargin(cell2mat(varargin)~=0))))./diag(hessianMatrix));
        end
        % For DSL model
    case 'HDSL'
        % Number of parameters
        Np = length(cell2mat(varargin(cell2mat(varargin)~=0)));
        % Calculate standard deviation of the data (not needed)
        stDevData = sqrt(1/(length(x)-length(Np)) * sum((z-fitVals).^2));
        % degree of variation of parameter to calculate sensitivity
        del = 0.001;
        % Generate and solve global optimisation problem for confidence regions
        % based on isotherm model
        qs1 = varargin{1};
        qs2 = varargin{2};
        b01 = varargin{3};
        b02 = varargin{4};
        delU1 = varargin{5};
        delU2 = varargin{6};
        qsH = varargin{7};
        b0H = varargin{8};
        delUH = varargin{9};
        switch fittingMethod
            case 'WSS'
                % Create empty sensitivity matrix
                sensitivityMatrix = zeros(length(x),9);
                % Calculate sensitivity at every data point for each parameter
                for jj = 1:9
                    for kk = 1:length(x)
                        if jj == 1
                            sensitivityMatrix(kk,jj) = (((((1+del)*qs1)*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + qs2*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))))) ...
                                + qsH*(b0H*x(kk)*exp(delUH/(8.314*y(kk))))) - fitVals(kk))/(del*qs1);
                        elseif jj == 2
                            sensitivityMatrix(kk,jj) = (((qs1*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + ((1+del)*qs2)*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))))) ...
                                + qsH*(b0H*x(kk)*exp(delUH/(8.314*y(kk))))) - fitVals(kk))/(del*qs2);
                        elseif jj == 3
                            sensitivityMatrix(kk,jj) = (((qs1*(((1+del)*b01)*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(((1+del)*b01)*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + qs2*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))))) ...
                                + qsH*(b0H*x(kk)*exp(delUH/(8.314*y(kk))))) - fitVals(kk))/(del*b01);
                        elseif jj == 4
                            sensitivityMatrix(kk,jj) = (((qs1*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + qs2*(((1+del)*b02)*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(((1+del)*b02)*x(kk)*exp(delU2/(8.314*y(kk)))))) ...
                                + qsH*(b0H*x(kk)*exp(delUH/(8.314*y(kk))))) - fitVals(kk))/(del*b02);
                        elseif jj == 5
                            sensitivityMatrix(kk,jj) = ((((qs1)*(b01*x(kk)*exp((1+del)*delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp((1+del)*delU1/(8.314*y(kk))))) ...
                                + qs2*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))))) ...
                                + qsH*(b0H*x(kk)*exp(delUH/(8.314*y(kk))))) - fitVals(kk))/(del*delU1);
                        elseif jj == 6
                            sensitivityMatrix(kk,jj) = ((((qs1)*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + qs2*(b02*x(kk)*exp((1+del)*delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp((1+del)*delU2/(8.314*y(kk)))))) ...
                                + qsH*(b0H*x(kk)*exp(delUH/(8.314*y(kk))))) - fitVals(kk))/(del*delU2);
                        elseif jj == 7
                            sensitivityMatrix(kk,jj) = ((((qs1)*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + qs2*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))))) ...
                                + (1+del)*qsH*(b0H*x(kk)*exp(delUH/(8.314*y(kk))))) - fitVals(kk))/(del*qsH);
                        elseif jj == 8
                            sensitivityMatrix(kk,jj) = (((qs1*((b01)*x(kk)*exp(delU1/(8.314*y(kk))))/(1+((b01)*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + qs2*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))))) ...
                                + qsH*((1+del)*b0H*x(kk)*exp(delUH/(8.314*y(kk))))) - fitVals(kk))/(del*b0H);
                        else
                            sensitivityMatrix(kk,jj) = ((((qs1)*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + qs2*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))))) ...
                                + qsH*(b0H*x(kk)*exp((1+del)*delUH/(8.314*y(kk))))) - fitVals(kk))/(del*delUH);
                        end
                    end
                end
                % estimated Hessian Matrix for the data set (Non-linear parameter estimation
                % by Yonathan Bard (1974) pg. 178
                hessianMatrix = 1/stDevData^2*transpose(sensitivityMatrix)*sensitivityMatrix;
                % Confidence range given by chi squared distribution at Np degrees
                % of freedom (independent parameter conf intervals)
                conRange95 = sqrt(chi2inv(0.95,Np)./diag(hessianMatrix));
            case 'MLE'
                Np = 9;
                Nt = length(x);
                del = 1e-6;
                dlogMLE = [];
                d2logMLE = [];
                deltaplus1mat = eye(Np).*(del);
                deltamat = eye(Np).*(del);
                partemp = [qs1./isoRef(1), qs2./isoRef(2),b01./isoRef(3),b02./isoRef(4),delU1./isoRef(5),delU2./isoRef(6),...
                    qsH./isoRef(7),b0H./isoRef(8),delUH./isoRef(9)];
                logMLE = @(par) -generateMLEfun(x, y, z, 1, 'HDSL', isoRef, par(1), par(2), par(3), ...
                    par(4), par(5), par(6), par(7), par(8), par(9));
                
                for jj = 1:Np
                    for kk = 1:Np
                        partempnumj = partemp.*(1+deltaplus1mat(jj,:));
                        partempdenj = partemp.*deltamat(jj,:);
                        partempnumk = partemp.*(1+deltaplus1mat(kk,:));
                        partempdenk = partemp.*deltamat(kk,:);
                        partempnumjk = partemp.*(1+deltaplus1mat(jj,:) + deltaplus1mat(kk,:));
                        d2logMLE(jj,kk) = ((logMLE(partempnumjk)-logMLE(partempnumk))-(logMLE(partempnumj)-logMLE(partemp)))./(partempdenj(jj).*isoRef(jj).*partempdenk(kk).*isoRef(kk));
                    end
                end
                hessianMatrix =  -d2logMLE;
                conRange95 = sqrt(chi2inv(0.95,length(cell2mat(varargin(cell2mat(varargin)~=0))))./diag(hessianMatrix));
        end
    case 'TSL'
        % Number of parameters
        Np = length(cell2mat(varargin(cell2mat(varargin)~=0)));
        % Calculate standard deviation of the data (not needed)
        stDevData = sqrt(1/(length(x)-length(Np)) * sum((z-fitVals).^2));
        % degree of variation of parameter to calculate sensitivity
        del = 0.001;
        % Generate and solve global optimisation problem for confidence regions
        % based on isotherm model
        qs1 = varargin{1};
        qs2 = varargin{2};
        b01 = varargin{3};
        b02 = varargin{4};
        delU1 = varargin{5};
        delU2 = varargin{6};
        qs3 = varargin{7};
        b03 = varargin{8};
        delU3 = varargin{9};
        switch fittingMethod
            case 'WSS'
                % Create empty sensitivity matrix
                sensitivityMatrix = zeros(length(x),9);
                % Calculate sensitivity at every data point for each parameter
                for jj = 1:9
                    for kk = 1:length(x)
                        if jj == 1
                            sensitivityMatrix(kk,jj) = (((((1+del)*qs1)*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + qs2*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))))) ...
                                + qs3*(b03*x(kk)*exp(delU3/(8.314*y(kk))))./(1+(b03*x(kk)*exp(delU3/(8.314*y(kk)))))) - fitVals(kk))/(del*qs1);
                        elseif jj == 2
                            sensitivityMatrix(kk,jj) = (((qs1*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + ((1+del)*qs2)*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))))) ...
                                + qs3*(b03*x(kk)*exp(delU3/(8.314*y(kk))))./(1+(b03*x(kk)*exp(delU3/(8.314*y(kk)))))) - fitVals(kk))/(del*qs2);
                        elseif jj == 3
                            sensitivityMatrix(kk,jj) = (((qs1*(((1+del)*b01)*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(((1+del)*b01)*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + qs2*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))))) ...
                                + qs3*(b03*x(kk)*exp(delU3/(8.314*y(kk))))./(1+(b03*x(kk)*exp(delU3/(8.314*y(kk)))))) - fitVals(kk))/(del*b01);
                        elseif jj == 4
                            sensitivityMatrix(kk,jj) = (((qs1*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + qs2*(((1+del)*b02)*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(((1+del)*b02)*x(kk)*exp(delU2/(8.314*y(kk)))))) ...
                                + qs3*(b03*x(kk)*exp(delU3/(8.314*y(kk))))./(1+(b03*x(kk)*exp(delU3/(8.314*y(kk)))))) - fitVals(kk))/(del*b02);
                        elseif jj == 5
                            sensitivityMatrix(kk,jj) = ((((qs1)*(b01*x(kk)*exp((1+del)*delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp((1+del)*delU1/(8.314*y(kk))))) ...
                                + qs2*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))))) ...
                                + qs3*(b03*x(kk)*exp(delU3/(8.314*y(kk))))./(1+(b03*x(kk)*exp(delU3/(8.314*y(kk)))))) - fitVals(kk))/(del*delU1);
                        elseif jj == 6
                            sensitivityMatrix(kk,jj) = ((((qs1)*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + qs2*(b02*x(kk)*exp((1+del)*delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp((1+del)*delU2/(8.314*y(kk)))))) ...
                                + qs3*(b03*x(kk)*exp(delU3/(8.314*y(kk))))./(1+(b03*x(kk)*exp(delU3/(8.314*y(kk)))))) - fitVals(kk))/(del*delU2);
                        elseif jj == 7
                            sensitivityMatrix(kk,jj) = ((((qs1)*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + qs2*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))))) ...
                                + (1+del)*qs3*(b03*x(kk)*exp(delU3/(8.314*y(kk))))./(1+(b03*x(kk)*exp(delU3/(8.314*y(kk)))))) - fitVals(kk))/(del*qs3);
                        elseif jj == 8
                            sensitivityMatrix(kk,jj) = (((qs1*((b01)*x(kk)*exp(delU1/(8.314*y(kk))))/(1+((b01)*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + qs2*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))))) ...
                                + qs3*((1+del)*b03*x(kk)*exp(delU3/(8.314*y(kk))))./(1+((1+del)*b03*x(kk)*exp(delU3/(8.314*y(kk)))))) - fitVals(kk))/(del*b03);
                        else
                            sensitivityMatrix(kk,jj) = ((((qs1)*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + qs2*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))))) ...
                                + qs3*(b03*x(kk)*exp((1+del)*delU3/(8.314*y(kk))))./(1+(b03*x(kk)*exp((1+del)*delU3/(8.314*y(kk)))))) - fitVals(kk))/(del*delU3);
                        end
                    end
                end
                % estimated Hessian Matrix for the data set (Non-linear parameter estimation
                % by Yonathan Bard (1974) pg. 178
                hessianMatrix = 1/stDevData^2*transpose(sensitivityMatrix)*sensitivityMatrix;
                % Confidence range given by chi squared distribution at Np degrees
                % of freedom (independent parameter conf intervals)
                conRange95 = sqrt(chi2inv(0.95,Np)./diag(hessianMatrix));
            case 'MLE'
                Nt = length(x);
                del = 1e-6;
                dlogMLE = [];
                d2logMLE = [];
                deltaplus1mat = eye(Np).*(del);
                deltamat = eye(Np).*(del);
                partemp = [qs1./isoRef(1), qs2./isoRef(2),b01./isoRef(3),b02./isoRef(4),delU1./isoRef(5),delU2./isoRef(6),...
                    qs3./isoRef(7),b03./isoRef(8),delU3./isoRef(9)];
                logMLE = @(par) -generateMLEfun(x, y, z, 1, 'TSL', isoRef, par(1), par(2), par(3), ...
                    par(4), par(5), par(6), par(7), par(8), par(9));
                
                for jj = 1:Np
                    for kk = 1:Np
                        partempnumj = partemp.*(1+deltaplus1mat(jj,:));
                        partempdenj = partemp.*deltamat(jj,:);
                        partempnumk = partemp.*(1+deltaplus1mat(kk,:));
                        partempdenk = partemp.*deltamat(kk,:);
                        partempnumjk = partemp.*(1+deltaplus1mat(jj,:) + deltaplus1mat(kk,:));
                        d2logMLE(jj,kk) = ((logMLE(partempnumjk)-logMLE(partempnumk))-(logMLE(partempnumj)-logMLE(partemp)))./(partempdenj(jj).*isoRef(jj).*partempdenk(kk).*isoRef(kk));
                    end
                end
                hessianMatrix =  -d2logMLE;
                conRange95 = sqrt(chi2inv(0.95,Np)./diag(hessianMatrix));
        end
        % for DSS model
    case 'DSS'
        % Number of parameters
        Np = length(cell2mat(varargin(cell2mat(varargin)~=0)));
        % Calculate standard deviation of the data (not needed)
        stDevData = sqrt(1/(length(x)-length(Np)) * sum((z-fitVals).^2));
        % degree of variation of parameter to calculate sensitivity
        del = 0.000001;
        % Generate and solve global optimisation problem for confidence regions
        % based on isotherm model
        qs1 = varargin{1};
        qs2 = varargin{2};
        b01 = varargin{3};
        b02 = varargin{4};
        delU1 = varargin{5};
        delU2 = varargin{6};
        gamma = varargin{7};
        switch fittingMethod
            case 'WSS'
                % Create empty sensitivity matrix
                sensitivityMatrix = zeros(length(x),7);
                % Calculate sensitivity at every data point for each parameter
                for jj = 1:7
                    for kk = 1:length(x)
                        if jj == 1
                            sensitivityMatrix(kk,jj) = ((((1+del)*qs1)*(b01*x(kk)*exp(delU1/(8.314*y(kk))))^gamma/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^gamma) ...
                                + qs2*(b02*x(kk)*exp(delU2/(8.314*y(kk))))^gamma/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk))))^gamma)) - fitVals(kk))/(del*qs1);
                        elseif jj == 2
                            sensitivityMatrix(kk,jj) = ((qs1*(b01*x(kk)*exp(delU1/(8.314*y(kk))))^gamma/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^gamma) ...
                                + ((1+del)*qs2)*(b02*x(kk)*exp(delU2/(8.314*y(kk))))^gamma/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk))))^gamma)) - fitVals(kk))/(del*qs2);
                        elseif jj == 3
                            sensitivityMatrix(kk,jj) = ((qs1*(((1+del)*b01)*x(kk)*exp(delU1/(8.314*y(kk))))^gamma/(1+(((1+del)*b01)*x(kk)*exp(delU1/(8.314*y(kk))))^gamma) ...
                                + qs2*(b02*x(kk)*exp(delU2/(8.314*y(kk))))^gamma/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk))))^gamma)) - fitVals(kk))/(del*b01);
                        elseif jj == 4
                            sensitivityMatrix(kk,jj) = ((qs1*(b01*x(kk)*exp(delU1/(8.314*y(kk))))^gamma/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^gamma) ...
                                + qs2*(((1+del)*b02)*x(kk)*exp(delU2/(8.314*y(kk))))^gamma/(1+(((1+del)*b02)*x(kk)*exp(delU2/(8.314*y(kk))))^gamma)) - fitVals(kk))/(del*b02);
                        elseif jj == 5
                            sensitivityMatrix(kk,jj) = ((qs1*(b01*x(kk)*exp((1+del)*delU1/(8.314*y(kk))))^gamma/(1+(b01*x(kk)*exp((1+del)*delU1/(8.314*y(kk))))^gamma) ...
                                + qs2*((b02)*x(kk)*exp(delU2/(8.314*y(kk))))^gamma/(1+((b02)*x(kk)*exp(delU2/(8.314*y(kk))))^gamma)) - fitVals(kk))/(del*delU1);
                        elseif jj == 6
                            sensitivityMatrix(kk,jj) = ((qs1*(b01*x(kk)*exp(delU1/(8.314*y(kk))))^gamma/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^gamma) ...
                                + qs2*((b02)*x(kk)*exp((1+del)*delU2/(8.314*y(kk))))^gamma/(1+((b02)*x(kk)*exp((1+del)*delU2/(8.314*y(kk))))^gamma)) - fitVals(kk))/(del*delU2);
                        else
                            sensitivityMatrix(kk,jj) = ((qs1*(b01*x(kk)*exp(delU1/(8.314*y(kk))))^((1+del)*gamma)/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^((1+del)*gamma)) ...
                                + qs2*((b02)*x(kk)*exp(delU2/(8.314*y(kk))))^((1+del)*gamma)/(1+((b02)*x(kk)*exp(delU2/(8.314*y(kk))))^((1+del)*gamma))) - fitVals(kk))/(del*gamma);
                        end
                    end
                end
                % Hessian Matrix for the data set (Non-linear parameter estimation
                % by Yonathan Bard (1974) pg. 178)
                hessianMatrix = 1/stDevData^2*transpose(sensitivityMatrix)*sensitivityMatrix;
                % Confidence range given by chi squared distribution at Np degrees
                % of freedom (independent parameter conf intervals)
                conRange95 = sqrt(chi2inv(0.95,Np)./diag(hessianMatrix));
            case 'MLE'
                Np = 7;
                Nt = length(x);
                del = 1e-6;
                dlogMLE = [];
                d2logMLE = [];
                deltaplus1mat = eye(Np).*(del);
                deltamat = eye(Np).*(del);
                partemp = [qs1./isoRef(1), qs2./isoRef(2),b01./isoRef(3),b02./isoRef(4),delU1./isoRef(5),delU2./isoRef(6),...
                    gamma./isoRef(7)];
                logMLE = @(par) -generateMLEfun(x, y, z, 1, 'DSS', isoRef, par(1), par(2), par(3), ...
                    par(4), par(5), par(6), par(7));
                
                for jj = 1:Np
                    for kk = 1:Np
                        partempnumj = partemp.*(1+deltaplus1mat(jj,:));
                        partempdenj = partemp.*deltamat(jj,:);
                        partempnumk = partemp.*(1+deltaplus1mat(kk,:));
                        partempdenk = partemp.*deltamat(kk,:);
                        partempnumjk = partemp.*(1+deltaplus1mat(jj,:) + deltaplus1mat(kk,:));
                        d2logMLE(jj,kk) = ((logMLE(partempnumjk)-logMLE(partempnumk))-(logMLE(partempnumj)-logMLE(partemp)))./(partempdenj(jj).*isoRef(jj).*partempdenk(kk).*isoRef(kk));
                    end
                end
                hessianMatrix =  -d2logMLE;
                conRange95 = sqrt(chi2inv(0.95,length(cell2mat(varargin(cell2mat(varargin)~=0))))./diag(hessianMatrix));
        end
        % for DSS model
    case 'TOTH'
        % Number of parameters
        Np = length(cell2mat(varargin(cell2mat(varargin)~=0)));
        % Calculate standard deviation of the data (not needed)
        stDevData = sqrt(1/(length(x)-length(Np)) * sum((z-fitVals).^2));
        % degree of variation of parameter to calculate sensitivity
        del = 0.000001;
        % Generate and solve global optimisation problem for confidence regions
        % based on isotherm model
        qs1 = varargin{1};
        qs2 = varargin{2};
        b01 = varargin{3};
        b02 = varargin{4};
        delU1 = varargin{5};
        delU2 = varargin{6};
        toth = varargin{7};
        switch fittingMethod
            case 'MLE'
                % Create empty sensitivity matrix
                sensitivityMatrix = zeros(length(x),7);
                % Calculate sensitivity at every data point for each parameter
                for jj = 1:7
                    for kk = 1:length(x)
                        if jj == 1
                            sensitivityMatrix(kk,jj) = ((((1+del)*qs1)*b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^toth).^(1./toth) ...
                                - fitVals(kk))/(del*qs1);
                        elseif jj == 2
                            sensitivityMatrix(kk,jj) = (((1+del)*qs1*b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^toth).^(1./toth) ...
                                - fitVals(kk))/(del*qs1);
                        elseif jj == 3
                            sensitivityMatrix(kk,jj) = ((qs1*(1+del)*b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^toth).^(1./toth) ...
                                - fitVals(kk))/(del*b01);
                        elseif jj == 4
                            sensitivityMatrix(kk,jj) = ((qs1*(1+del)*b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^toth).^(1./toth) ...
                                - fitVals(kk))/(del*b01);
                        elseif jj == 5
                            sensitivityMatrix(kk,jj) = ((qs1*b01*x(kk)*exp((1+del)*delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^toth).^(1./toth) ...
                                - fitVals(kk))/(del*delU1);
                        elseif jj == 6
                            sensitivityMatrix(kk,jj) = ((qs1*b01*x(kk)*exp((1+del)*delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^toth).^(1./toth) ...
                                - fitVals(kk))/(del*delU1);
                        else
                            sensitivityMatrix(kk,jj) = ((qs1*b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^(1+del)*toth).^(1./(1+del)*toth) ...
                                - fitVals(kk))/(del*toth);
                        end
                    end
                end
                % Hessian Matrix for the data set (Non-linear parameter estimation
                % by Yonathan Bard (1974) pg. 178)
                hessianMatrix = 1/stDevData^2*transpose(sensitivityMatrix)*sensitivityMatrix;
                % Confidence range given by chi squared distribution at Np degrees
                % of freedom (independent parameter conf intervals)
                conRange95 = sqrt(chi2inv(0.95,Np)./diag(hessianMatrix));
            case 'WSS'
                Np = 7;
                Nt = length(x);
                del = 1e-6;
                dlogMLE = [];
                d2logMLE = [];
                deltaplus1mat = eye(Np).*(del);
                deltamat = eye(Np).*(del);
                partemp = [qs1./isoRef(1), b01./isoRef(3), delU1./isoRef(5), toth./isoRef(7)];
                logMLE = @(par) -generateMLEfun(x, y, z, 1, 'TOTH', isoRef, par(1), 0, par(2), ...
                    0, par(3), 0, par(4));
                
                for jj = 1:Np
                    for kk = 1:Np
                        partempnumj = partemp.*(1+deltaplus1mat(jj,:));
                        partempdenj = partemp.*deltamat(jj,:);
                        partempnumk = partemp.*(1+deltaplus1mat(kk,:));
                        partempdenk = partemp.*deltamat(kk,:);
                        partempnumjk = partemp.*(1+deltaplus1mat(jj,:) + deltaplus1mat(kk,:));
                        d2logMLE(jj,kk) = ((logMLE(partempnumjk)-logMLE(partempnumk))-(logMLE(partempnumj)-logMLE(partemp)))./(partempdenj(jj).*isoRef(jj).*partempdenk(kk).*isoRef(kk));
                    end
                end
                hessianMatrix =  -d2logMLE;
                conRange95 = sqrt(chi2inv(0.95, length(cell2mat(varargin(cell2mat(varargin)~=0))))./diag(hessianMatrix));
        end
        % for Virial model
    case 'VIRIAL'
        % Number of parameters
        Np = length(cell2mat(varargin(cell2mat(varargin)~=0)));
        % Calculate standard deviation of the data (not needed)
        stDevData = sqrt(1/(length(x)-length(Np)) * sum((log(x)-fitVals').^2));
        % degree of variation of parameter to calculate sensitivity
        del = 0.000001;
        % Generate and solve global optimisation problem for confidence regions
        % based on isotherm model
        a0 = varargin{1};
        a1 = varargin{2};
        a2 = varargin{3};
        a3 = varargin{4};
        b0 = varargin{5};
        b1 = varargin{6};
        b2 = varargin{7};
        b3 = varargin{8};
        
        parameters = [a0, a1, a2, a3, b0, b1, b2, b3];
        switch fittingMethod
            case 'WSS'
                % Create empty sensitivity matrix
                sensitivityMatrix = zeros(length(x),8);
                % Calculate sensitivity at every data point for each parameter
                for jj = 1:8
                    for kk = 1:length(x)
                        if jj == 1
                            sensitivityMatrix(kk,jj) = (log(z(kk)) + 1/y(kk).*((1+del)*parameters(1) + parameters(2)*z(kk) + ...
                                parameters(3)*z(kk)^2 + parameters(4)*z(kk)^3)  + parameters(5) ...
                                + parameters(6)*z(kk)+ parameters(7)*z(kk).^2+ parameters(8)*z(kk).^3 - fitVals(kk))/(del*parameters(1));
                        elseif jj == 2
                            sensitivityMatrix(kk,jj) = (log(z(kk)) + 1/y(kk).*(parameters(1) + (1+del)*parameters(2)*z(kk) + ...
                                parameters(3)*z(kk)^2 + parameters(4)*z(kk)^3)  + parameters(5) ...
                                + parameters(6)*z(kk)+ parameters(7)*z(kk).^2+ parameters(8)*z(kk).^3 - fitVals(kk))/(del*parameters(2));
                        elseif jj == 3
                            sensitivityMatrix(kk,jj) = (log(z(kk)) + 1/y(kk).*(parameters(1) + parameters(2)*z(kk) + ...
                                (1+del)*parameters(3)*z(kk)^2 + parameters(4)*z(kk)^3)  + parameters(5) ...
                                + parameters(6)*z(kk)+ parameters(7)*z(kk).^2+ parameters(8)*z(kk).^3 - fitVals(kk))/(del*parameters(3));
                        elseif jj == 4
                            sensitivityMatrix(kk,jj) = (log(z(kk)) + 1/y(kk).*(parameters(1) + parameters(2)*z(kk) + ...
                                parameters(3)*z(kk)^2 + (1+del)*parameters(4)*z(kk)^3)  + parameters(5) ...
                                + parameters(6)*z(kk)+ parameters(7)*z(kk).^2+ parameters(8)*z(kk).^3 - fitVals(kk))/(del*parameters(4));
                        elseif jj == 5
                            sensitivityMatrix(kk,jj) = (log(z(kk)) + 1/y(kk).*(parameters(1) + parameters(2)*z(kk) + ...
                                parameters(3)*z(kk)^2 + parameters(4)*z(kk)^3)  + (1+del)*parameters(5) ...
                                + parameters(6)*z(kk)+ parameters(7)*z(kk).^2+ parameters(8)*z(kk).^3 - fitVals(kk))/(del*parameters(5));
                        elseif jj == 6
                            sensitivityMatrix(kk,jj) = (log(z(kk)) + 1/y(kk).*(parameters(1) + parameters(2)*z(kk) + ...
                                parameters(3)*z(kk)^2 + parameters(4)*z(kk)^3)  + parameters(5) ...
                                + (1+del)*parameters(6)*z(kk)+ parameters(7)*z(kk).^2+ parameters(8)*z(kk).^3 - fitVals(kk))/(del*parameters(6));
                        elseif jj == 7
                            sensitivityMatrix(kk,jj) = (log(z(kk)) + 1/y(kk).*(parameters(1) + parameters(2)*z(kk) + ...
                                parameters(3)*z(kk)^2 + parameters(4)*z(kk)^3)  + parameters(5) ...
                                + parameters(6)*z(kk)+ (1+del)*parameters(7)*z(kk).^2+ parameters(8)*z(kk).^3 - fitVals(kk))/(del*parameters(7));
                        else
                            sensitivityMatrix(kk,jj) = (log(z(kk)) + 1/y(kk).*(parameters(1) + parameters(2)*z(kk) + ...
                                parameters(3)*z(kk)^2 + parameters(4)*z(kk)^3)  + parameters(5) ...
                                + parameters(6)*z(kk)+ parameters(7)*z(kk).^2+ (1+del)*parameters(8)*z(kk).^3 - fitVals(kk))/(del*parameters(8));
                        end
                    end
                end
                % Hessian Matrix for the data set (Non-linear parameter estimation
                % by Yonathan Bard (1974) pg. 178)
                hessianMatrix = 1/stDevData^2*transpose(sensitivityMatrix)*sensitivityMatrix;
                % Confidence range given by chi squared distribution at Np degrees
                % of freedom (independent parameter conf intervals)
                conRange95 = sqrt(chi2inv(0.95,Np)./diag(hessianMatrix));
            case 'MLE'
                Nt = length(x);
                del = 1e-6;
                dlogMLE = [];
                d2logMLE = [];
                deltaplus1mat = eye(Np).*(del);
                deltamat = eye(Np).*(del);
                partemp = [a0./isoRef(1), a1./isoRef(2),a2./isoRef(3),a3./isoRef(4),b0./isoRef(5),b1./isoRef(6)];
                logMLE = @(par) -generateMLEfun(x, y, z, 1, 'VIRIAL', isoRef, par(1), par(2), par(3), ...
                    par(4), par(5), par(6), 0, 0);
                
                for jj = 1:Np
                    for kk = 1:Np
                        partempnumj = partemp.*(1+deltaplus1mat(jj,:));
                        partempdenj = partemp.*deltamat(jj,:);
                        partempnumk = partemp.*(1+deltaplus1mat(kk,:));
                        partempdenk = partemp.*deltamat(kk,:);
                        partempnumjk = partemp.*(1+deltaplus1mat(jj,:) + deltaplus1mat(kk,:));
                        d2logMLE(jj,kk) = ((logMLE(partempnumjk)-logMLE(partempnumk))-(logMLE(partempnumj)-logMLE(partemp)))./(partempdenj(jj).*isoRef(jj).*partempdenk(kk).*isoRef(kk));
                    end
                end
                hessianMatrix =  -d2logMLE;
                conRange95 = sqrt(chi2inv(0.95,Np)./diag(hessianMatrix));
        end
    case 'VIRIAL2'
        % Number of parameters
        Np = length(cell2mat(varargin));
        % Calculate standard deviation of the data (not needed)
        stDevData = sqrt(1/(length(x)-length(Np)) * sum((log(x)-fitVals').^2));
        % degree of variation of parameter to calculate sensitivity
        del = 0.000001;
        % Generate and solve global optimisation problem for confidence regions
        % based on isotherm model
        a0 = varargin{1};
        a1 = varargin{2};
        a2 = varargin{3};
        a3 = varargin{4};
        b0 = varargin{5};
        b1 = varargin{6};
        b2 = varargin{7};
        b3 = varargin{8};
        
        parameters = [a0, a1, a2, a3, b0, b1, b2, b3];
        switch fittingMethod
            case 'WSS'
                % Create empty sensitivity matrix
                sensitivityMatrix = zeros(length(x),8);
                % Calculate sensitivity at every data point for each parameter
                for jj = 1:8
                    for kk = 1:length(x)
                        if jj == 1
                            sensitivityMatrix(kk,jj) = (log(z(kk)) + 1/y(kk).*((1+del)*parameters(1) + parameters(2)*z(kk) + ...
                                parameters(3)*z(kk)^2 + parameters(4)*z(kk)^3)  + parameters(5) ...
                                + parameters(6)*z(kk)+ parameters(7)*z(kk).^2+ parameters(8)*z(kk).^3 - fitVals(kk))/(del*parameters(1));
                        elseif jj == 2
                            sensitivityMatrix(kk,jj) = (log(z(kk)) + 1/y(kk).*(parameters(1) + (1+del)*parameters(2)*z(kk) + ...
                                parameters(3)*z(kk)^2 + parameters(4)*z(kk)^3)  + parameters(5) ...
                                + parameters(6)*z(kk)+ parameters(7)*z(kk).^2+ parameters(8)*z(kk).^3 - fitVals(kk))/(del*parameters(2));
                        elseif jj == 3
                            sensitivityMatrix(kk,jj) = (log(z(kk)) + 1/y(kk).*(parameters(1) + parameters(2)*z(kk) + ...
                                (1+del)*parameters(3)*z(kk)^2 + parameters(4)*z(kk)^3)  + parameters(5) ...
                                + parameters(6)*z(kk)+ parameters(7)*z(kk).^2+ parameters(8)*z(kk).^3 - fitVals(kk))/(del*parameters(3));
                        elseif jj == 4
                            sensitivityMatrix(kk,jj) = (log(z(kk)) + 1/y(kk).*(parameters(1) + parameters(2)*z(kk) + ...
                                parameters(3)*z(kk)^2 + (1+del)*parameters(4)*z(kk)^3)  + parameters(5) ...
                                + parameters(6)*z(kk)+ parameters(7)*z(kk).^2+ parameters(8)*z(kk).^3 - fitVals(kk))/(del*parameters(4));
                        elseif jj == 5
                            sensitivityMatrix(kk,jj) = (log(z(kk)) + 1/y(kk).*(parameters(1) + parameters(2)*z(kk) + ...
                                parameters(3)*z(kk)^2 + parameters(4)*z(kk)^3)  + (1+del)*parameters(5) ...
                                + parameters(6)*z(kk)+ parameters(7)*z(kk).^2+ parameters(8)*z(kk).^3 - fitVals(kk))/(del*parameters(5));
                        elseif jj == 6
                            sensitivityMatrix(kk,jj) = (log(z(kk)) + 1/y(kk).*(parameters(1) + parameters(2)*z(kk) + ...
                                parameters(3)*z(kk)^2 + parameters(4)*z(kk)^3)  + parameters(5) ...
                                + (1+del)*parameters(6)*z(kk)+ parameters(7)*z(kk).^2+ parameters(8)*z(kk).^3 - fitVals(kk))/(del*parameters(6));
                        elseif jj == 7
                            sensitivityMatrix(kk,jj) = (log(z(kk)) + 1/y(kk).*(parameters(1) + parameters(2)*z(kk) + ...
                                parameters(3)*z(kk)^2 + parameters(4)*z(kk)^3)  + parameters(5) ...
                                + parameters(6)*z(kk)+ (1+del)*parameters(7)*z(kk).^2+ parameters(8)*z(kk).^3 - fitVals(kk))/(del*parameters(7));
                        else
                            sensitivityMatrix(kk,jj) = (log(z(kk)) + 1/y(kk).*(parameters(1) + parameters(2)*z(kk) + ...
                                parameters(3)*z(kk)^2 + parameters(4)*z(kk)^3)  + parameters(5) ...
                                + parameters(6)*z(kk)+ parameters(7)*z(kk).^2+ (1+del)*parameters(8)*z(kk).^3 - fitVals(kk))/(del*parameters(8));
                        end
                    end
                end
                % Hessian Matrix for the data set (Non-linear parameter estimation
                % by Yonathan Bard (1974) pg. 178)
                hessianMatrix = 1/stDevData^2*transpose(sensitivityMatrix)*sensitivityMatrix;
                % Confidence range given by chi squared distribution at Np degrees
                % of freedom (independent parameter conf intervals)
                conRange95 = sqrt(chi2inv(0.95,Np)./diag(hessianMatrix));
            case 'MLE'
                Nt = length(x);
                del = 1e-6;
                dlogMLE = [];
                d2logMLE = [];
                deltaplus1mat = eye(Np).*(del);
                deltamat = eye(Np).*(del);
                partemp = [a0./isoRef(1), a1./isoRef(2),a2./isoRef(3),a3./isoRef(4),b0./isoRef(5),b1./isoRef(6),b2./isoRef(7),b3./isoRef(8)];
                logMLE = @(par) -generateMLEfun(x, y, z, 1, 'VIRIAL2', isoRef, par(1), par(2), par(3), ...
                    par(4), par(5), par(6), par(7), par(8));
                
                for jj = 1:Np
                    for kk = 1:Np
                        partempnumj = partemp.*(1+deltaplus1mat(jj,:));
                        partempdenj = partemp.*deltamat(jj,:);
                        partempnumk = partemp.*(1+deltaplus1mat(kk,:));
                        partempdenk = partemp.*deltamat(kk,:);
                        partempnumjk = partemp.*(1+deltaplus1mat(jj,:) + deltaplus1mat(kk,:));
                        d2logMLE(jj,kk) = ((logMLE(partempnumjk)-logMLE(partempnumk))-(logMLE(partempnumj)-logMLE(partemp)))./(partempdenj(jj).*isoRef(jj).*partempdenk(kk).*isoRef(kk));
                    end
                end
                hessianMatrix =  -d2logMLE;
                conRange95 = sqrt(chi2inv(0.95,Np)./diag(hessianMatrix));
        end
end
end