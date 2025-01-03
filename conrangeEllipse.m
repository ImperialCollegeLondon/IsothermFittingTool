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
% degree of variation of parameter to calculate sensitivity
del = 1e-5;
switch isothermModel
    % Calculate error for Statistical isotherm model for zeolites
    case 'STATZ'
        % Calculate standard deviation of the data (not needed)
        Np =4;
        stDevData = sqrt(1/(length(x)-Np) * sum((z-fitVals).^2));
        omega = varargin{1};
        beta = varargin{2};
        b01 = varargin{3};
        delU1 = varargin{4};
        vc = varargin{5};
        vm = varargin{6};
        Nt = length(x);

        dlogMLE = [];
        d2logMLE = [];
        deltaplus1mat = eye(Np).*(del);
        deltamat = eye(Np).*(del);
        partemp = [omega./isoRef(1), beta./isoRef(2),b01./isoRef(3),delU1./isoRef(4)];
        logMLE = @(par) -generateMLEfun(x, y, z, 1, 'STATZ', isoRef, par(1), par(2), par(3), par(4), vc, vm);

        for jj = 2:Np
            for kk = 2:Np
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
        conRange95 = sqrt(chi2inv(0.95,Np)./diag(hessianMatrix));
    case 'STATZE'
        % Calculate standard deviation of the data (not needed)
        Np = 5;
        stDevData = sqrt(1/(length(x)-Np) * sum((z-fitVals).^2));
        omega = varargin{1};
        beta = varargin{2};
        b01 = varargin{3};
        delU1 = varargin{4};
        ebyk = varargin{5};
        vc = varargin{6};
        vm = varargin{7};
        Nt = length(x);

        dlogMLE = [];
        d2logMLE = [];
        deltaplus1mat = eye(Np).*(del);
        deltamat = eye(Np).*(del);
        partemp = [omega./isoRef(1), beta./isoRef(2),b01./isoRef(3),delU1./isoRef(4),ebyk./isoRef(5)];
        logMLE = @(par) -generateMLEfun(x, y, z, 1, 'STATZE', isoRef, par(1), par(2), par(3), par(4), par(5), vc, vm);

        for jj = 2:Np
            for kk = 2:Np
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
        conRange95 = sqrt(chi2inv(0.95,Np)./diag(hessianMatrix));
    case 'STATZGO'
        % Calculate standard deviation of the data (not needed)
        Np =7;
        stDevData = sqrt(1/(length(x)-Np) * sum((z-fitVals).^2));
        omega = varargin{1};
        beta = varargin{2};
        b01 = varargin{3};
        delU1 = varargin{4};
        delU2 = varargin{5};
        kgate = varargin{6};
        cgate = varargin{7};
        vc = varargin{8};
        vm = varargin{9};
        Nt = length(x);

        dlogMLE = [];
        d2logMLE = [];
        deltaplus1mat = eye(Np).*(del);
        deltamat = eye(Np).*(del);
        partemp = [omega./isoRef(1), beta./isoRef(2),b01./isoRef(3),delU1./isoRef(4),delU2./isoRef(5),kgate./isoRef(6),cgate./isoRef(7)];
        logMLE = @(par) -generateMLEfun(x, y, z, 1, 'STATZGO', isoRef, par(1), par(2), par(3), par(4),par(5), par(6), par(7), vc, vm);

        for jj = 1:Np
            for kk = 2:Np
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
        hessianMatrix =  abs(d2logMLE);
        conRange95 = sqrt(chi2inv(0.95,Np)./diag(hessianMatrix));
    case 'SSLSTA'
        % Calculate standard deviation of the data (not needed)
        Np = 9;
        stDevData = sqrt(1/(length(x)-Np) * sum((z-fitVals).^2));
        qsNP = varargin{1};
        b01NP = varargin{2};
        delU1NP = varargin{3};
        qsLP = varargin{4};
        b01LP = varargin{5};
        delU1LP = varargin{6};
        kgate = varargin{7};
        cgate = varargin{8};
        sval = varargin{9};

        Nt = length(x);

        dlogMLE = [];
        d2logMLE = [];
        deltaplus1mat = eye(Np).*(del);
        deltamat = eye(Np).*(del);
        partemp = [qsNP./isoRef(1), b01NP./isoRef(2),delU1NP./isoRef(3),qsLP./isoRef(4),b01LP./isoRef(5),delU1LP./isoRef(6),kgate./isoRef(7),cgate./isoRef(8),sval./isoRef(9)];
        logMLE = @(par) -generateMLEfun(x, y, z, 1, 'SSLSTA', isoRef, par(1), par(2), par(3), par(4),par(5), par(6), par(7), par(8), par(9));

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
        hessianMatrix =  abs(d2logMLE);
        conRange95 = sqrt(chi2inv(0.95,Np)./diag(hessianMatrix));
    case 'STATSTA2'
        % Calculate standard deviation of the data (not needed)
        Np = 9;
        stDevData = sqrt(1/(length(x)-Np) * sum((z-fitVals).^2));
        omega = varargin{1};
        beta = varargin{2};
        b01 = varargin{3};
        b02 = varargin{4};
        delU1 = varargin{5};
        delU2 = varargin{6};
        kgate = varargin{7};
        cgate = varargin{8};
        delvc = varargin{9};
        vc = varargin{10};
        vm = varargin{11};
        vc_old = varargin{12};
        Nt = length(x);
        cutoffq = varargin{13};

        isoRef = [isoRef(1) isoRef(2) isoRef(3) isoRef(3) isoRef(4) isoRef(5) isoRef(6) isoRef(7) isoRef(8)];

        dlogMLE = [];
        d2logMLE = [];
        deltaplus1mat = eye(Np).*(del);
        deltamat = eye(Np).*(del);
        partemp = [omega./isoRef(1), beta./isoRef(2),b01./isoRef(3), b02./isoRef(3), delU1./isoRef(4),delU2./isoRef(5),kgate./isoRef(6),cgate./isoRef(7),delvc./isoRef(8)];
        logMLE = @(par) -generateMLEfun(x, y, z, 1, 'STATSTA2', isoRef, par(1), par(2), par(3), par(4),par(5), par(6), par(7), par(8), par(9), vc, vm, vc_old, cutoffq);

        for jj = 1:Np-1
            for kk = 1:Np-1
                partempnumj = partemp.*(1+deltaplus1mat(jj+1,:));
                partempdenj = partemp.*deltamat(jj+1,:);
                partempnumk = partemp.*(1+deltaplus1mat(kk+1,:));
                partempdenk = partemp.*deltamat(kk+1,:);
                partempnumjk = partemp.*(1+deltaplus1mat(jj+1,:) + deltaplus1mat(kk+1,:));
                % Compute second derivative of logL for jj and kk
                d2logMLE(jj,kk) = ((logMLE(partempnumjk)-logMLE(partempnumk))-(logMLE(partempnumj)-logMLE(partemp)))./(partempdenj(jj).*isoRef(jj).*partempdenk(kk).*isoRef(kk)+1);
            end
        end
        % estimated Hessian Matrix for the data set (Non-linear parameter estimation
        % by Yonathan Bard (1974) pg. 178 (Eqn 7-5-17)
        hessianMatrix =  abs(d2logMLE);
        conRange95 = sqrt(chi2inv(0.95,size(hessianMatrix(:,1)))./diag(hessianMatrix));
    case 'STATSTA2gamma'
        % Calculate standard deviation of the data (not needed)
        Np = 10;
        stDevData = sqrt(1/(length(x)-Np) * sum((z-fitVals).^2));
        omega = varargin{1};
        beta = varargin{2};
        b01 = varargin{3};
        b02 = varargin{4};
        delU1 = varargin{5};
        delU2 = varargin{6};
        kgate = varargin{7};
        cgate = varargin{8};
        delvc = varargin{9};
        gamma = varargin{10};
        vc = varargin{11};
        vm = varargin{12};
        vc_old = varargin{13};
        Nt = length(x);
        cutoffq = varargin{14};

        isoRef = [isoRef(1) isoRef(2) isoRef(3) isoRef(3) isoRef(4) isoRef(5) isoRef(6) isoRef(7) isoRef(8) isoRef(9)];

        dlogMLE = [];
        d2logMLE = [];
        deltaplus1mat = eye(Np).*(del);
        deltamat = eye(Np).*(del);
        partemp = [omega./isoRef(1), beta./isoRef(2),b01./isoRef(3), b02./isoRef(3), delU1./isoRef(4),delU2./isoRef(5),kgate./isoRef(6),cgate./isoRef(7),delvc./isoRef(8), gamma./isoRef(9)];
        logMLE = @(par) -generateMLEfun(x, y, z, 1, 'STATSTA2gamma', isoRef, par(1), par(2), par(3), par(4),par(5), par(6), par(7), par(8), par(9), par(10), vc, vm, vc_old, cutoffq);

        for jj = 1:Np-1
            for kk = 1:Np-1
                partempnumj = partemp.*(1+deltaplus1mat(jj+1,:));
                partempdenj = partemp.*deltamat(jj+1,:);
                partempnumk = partemp.*(1+deltaplus1mat(kk+1,:));
                partempdenk = partemp.*deltamat(kk+1,:);
                partempnumjk = partemp.*(1+deltaplus1mat(jj+1,:) + deltaplus1mat(kk+1,:));
                % Compute second derivative of logL for jj and kk
                d2logMLE(jj,kk) = ((logMLE(partempnumjk)-logMLE(partempnumk))-(logMLE(partempnumj)-logMLE(partemp)))./(partempdenj(jj).*isoRef(jj).*partempdenk(kk).*isoRef(kk)+1);
            end
        end
        % estimated Hessian Matrix for the data set (Non-linear parameter estimation
        % by Yonathan Bard (1974) pg. 178 (Eqn 7-5-17)
        hessianMatrix =  abs(d2logMLE);
        conRange95 = sqrt(chi2inv(0.95,size(hessianMatrix(:,1)))./diag(hessianMatrix));
    case 'STATZSips'
        % Calculate standard deviation of the data (not needed)
        Np =5;
        stDevData = sqrt(1/(length(x)-Np) * sum((z-fitVals).^2));
        omega = varargin{1};
        beta = varargin{2};
        b01 = varargin{3};
        delU1 = varargin{4};
        gamma = varargin{5};
        vc = varargin{6};
        vm = varargin{7};
        Nt = length(x);

        dlogMLE = [];
        d2logMLE = [];
        deltaplus1mat = eye(Np).*(del);
        deltamat = eye(Np).*(del);
        partemp = [omega./isoRef(1), beta./isoRef(2),b01./isoRef(3),delU1./isoRef(4),gamma./isoRef(5)];
        logMLE = @(par) -generateMLEfun(x, y, z, 1, 'STATZSips', isoRef, par(1), par(2), par(3), par(4), par(5), vc, vm);

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
        conRange95 = sqrt(chi2inv(0.95,Np)./diag(hessianMatrix));
    case 'GAB'
        % Calculate standard deviation of the data (not needed)
        Np =5;
        stDevData = sqrt(1/(length(x)-Np) * sum((z-fitVals).^2));
        qs1  = varargin{1};
        parC = varargin{2};
        parD = varargin{3};
        parF = varargin{4};
        parG = varargin{5};
        Nt = length(x);

        dlogMLE = [];
        d2logMLE = [];
        deltaplus1mat = eye(Np).*(del);
        deltamat = eye(Np).*(del);
        partemp = [qs1./isoRef(1), parC./isoRef(2),parD./isoRef(3),parF./isoRef(4),parG./isoRef(5)];
        logMLE = @(par) -generateMLEfun(x, y, z, 1, 'GAB', isoRef, par(1), par(2), par(3), par(4), par(5));

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
        conRange95 = sqrt(chi2inv(0.95,Np)./diag(hessianMatrix));
    case 'UNIV6'
        % Calculate standard deviation of the data (not needed)
        Np =13;
        stDevData = sqrt(1/(length(x)-Np).* sum((z-fitVals').^2));
        qs =  varargin{1};
        a1 =  varargin{2};
        a2 =  varargin{3};
        a3 =  varargin{4};
        e01 = varargin{5};
        e02 = varargin{6};
        e03 = varargin{7};
        e04 = varargin{8};
        m1 =  varargin{9};
        m2 = varargin{10};
        m3 = varargin{11};
        m4 = varargin{12};
        ps = varargin{13};

        %         Create empty sensitivity matrix
        sensitivityMatrix = zeros(length(x),Np);
        del = 1e-3;
        % Calculate sensitivity at every data point for each parameter
        for jj = 1:Np
            for kk = 1:length(x)
                if jj == 1
                    sensitivityMatrix(kk,jj) = (computeUNIV6Loading(x(kk),y(kk), (1+del).*qs, a1, a2, a3, e01, e02, e03, e04, m1, m2, m3, m4, ps) - fitVals(kk))/(del*qs);
                elseif jj == 2
                    sensitivityMatrix(kk,jj) = (computeUNIV6Loading(x(kk),y(kk), qs, (1+del).*a1, a2, a3, e01, e02, e03, e04, m1, m2, m3, m4, ps) - fitVals(kk))/(del*a1);
                elseif jj == 3
                    sensitivityMatrix(kk,jj) = (computeUNIV6Loading(x(kk),y(kk), qs, a1, (1+del).*a2, a3, e01, e02, e03, e04, m1, m2, m3, m4, ps) - fitVals(kk))/(del*a2);
                elseif jj == 4
                    sensitivityMatrix(kk,jj) = (computeUNIV6Loading(x(kk),y(kk), qs, a1, a2, (1+del).*a3, e01, e02, e03, e04, m1, m2, m3, m4, ps) - fitVals(kk))/(del*a3);
                elseif jj == 5
                    sensitivityMatrix(kk,jj) = (computeUNIV6Loading(x(kk),y(kk), qs, a1, a2, a3, (1+del).*e01, e02, e03, e04, m1, m2, m3, m4, ps) - fitVals(kk))/(del*e01);
                elseif jj == 6
                    sensitivityMatrix(kk,jj) = (computeUNIV6Loading(x(kk),y(kk), qs, a1, a2, a3, e01, (1+del).*e02, e03, e04, m1, m2, m3, m4, ps) - fitVals(kk))/(del*e02);
                elseif jj == 7
                    sensitivityMatrix(kk,jj) = (computeUNIV6Loading(x(kk),y(kk), qs, a1, a2, a3, e01, e02, (1+del).*e03, e04, m1, m2, m3, m4, ps) - fitVals(kk))/(del*e03);
                elseif jj == 8
                    sensitivityMatrix(kk,jj) = (computeUNIV6Loading(x(kk),y(kk), qs, a1, a2, a3, e01, e02, e03, (1+del).*e04, m1, m2, m3, m4, ps) - fitVals(kk))/(del*e04);
                elseif jj == 9
                    sensitivityMatrix(kk,jj) = (computeUNIV6Loading(x(kk),y(kk), qs, a1, a2, a3, e01, e02, e03, e04, (1+del).*m1, m2, m3, m4, ps) - fitVals(kk))/(del*m1);
                elseif jj == 10
                    sensitivityMatrix(kk,jj) = (computeUNIV6Loading(x(kk),y(kk), qs, a1, a2, a3, e01, e02, e03, e04, m1, (1+del).*m2, m3, m4, ps) - fitVals(kk))/(del*m2);
                elseif jj == 11
                    sensitivityMatrix(kk,jj) = (computeUNIV6Loading(x(kk),y(kk), qs, a1, a2, a3, e01, e02, e03, e04, m1, m2, (1+del).*m3, m4, ps) - fitVals(kk))/(del*m3);
                elseif jj == 12
                    sensitivityMatrix(kk,jj) = (computeUNIV6Loading(x(kk),y(kk), qs, a1, a2, a3, e01, e02, e03, e04, m1, m2, m3, (1+del).*m4, ps) - fitVals(kk))/(del*m4);
                else
                    sensitivityMatrix(kk,jj) = (computeUNIV6Loading(x(kk),y(kk), qs, a1, a2, a3, e01, e02, e03, e04, m1, m2, m3, m4, (1+del).*ps) - fitVals(kk))/(del*ps);
                end
            end
        end
        % estimated Hessian Matrix for the data set (Non-linear parameter estimation
        % by Yonathan Bard (1974) pg. 178
        hessianMatrix = 1/stDevData.^2*transpose(sensitivityMatrix)*sensitivityMatrix;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         dlogMLE = [];
        %         d2logMLE = [];
        %         deltaplus1mat = eye(Np).*(del);
        %         deltamat = eye(Np).*(del);
        %         partemp = [qs./isoRef(1), a1./isoRef(2), a2./isoRef(3), a3./isoRef(4), ...
        %             e01./isoRef(5), e02./isoRef(6), e03./isoRef(7), e04./isoRef(8), ...
        %             m1./isoRef(9), m2./isoRef(10), m3./isoRef(11), m4./isoRef(12) , ps./isoRef(13) ];
        %         logMLE = @(par) -generateMLEfun(x, y, z, 1, 'UNIV6', isoRef, par(1), par(2), ...
        %                         par(3), par(4), par(5), par(6), par(7), par(8), par(9), par(10), par(11), par(12), par(13));
        %
        %         for jj = 1:Np
        %             for kk = 1:Np
        %                 partempnumj = partemp.*(1+deltaplus1mat(jj,:));
        %                 partempdenj = partemp.*deltamat(jj,:);
        %                 partempnumk = partemp.*(1+deltaplus1mat(kk,:));
        %                 partempdenk = partemp.*deltamat(kk,:);
        %                 partempnumjk = partemp.*(1+deltaplus1mat(jj,:) + deltaplus1mat(kk,:));
        %                 % Compute second derivative of logL for jj and kk
        %                 d2logMLE(jj,kk) = ((logMLE(partempnumjk)-logMLE(partempnumk))-(logMLE(partempnumj)-logMLE(partemp)))./(partempdenj(jj).*isoRef(jj).*partempdenk(kk).*isoRef(kk));
        %             end
        %         end
        %         % estimated Hessian Matrix for the data set (Non-linear parameter estimation
        %         % by Yonathan Bard (1974) pg. 178 (Eqn 7-5-17)
        %         hessianMatrix =  -d2logMLE;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        conRange95 = sqrt(chi2inv(0.95,Np)./diag(hessianMatrix));

        % Calculate error for Statistical isotherm model for zeolites
    case 'STATZGATE'
        % Calculate standard deviation of the data (not needed)
        Np = 7;
        stDevData = sqrt(1/(length(x)-Np) * sum((z-fitVals).^2));
        omega = varargin{1};
        beta = varargin{2};
        b01 = varargin{3};
        delU1 = varargin{4};
        kgate = varargin{5};
        cgate = varargin{6};
        gamma = varargin{7};
        vc1 = varargin{8};
        vc2 = varargin{9};
        vm = varargin{10};

        Nt = length(x);

        dlogMLE = [];
        d2logMLE = [];
        deltaplus1mat = eye(Np).*(del);
        deltamat = eye(Np).*(del);
        partemp = [omega./isoRef(1), beta./isoRef(2),b01./isoRef(3),delU1./isoRef(4), kgate./isoRef(5), cgate./isoRef(6), gamma./isoRef(7)];
        logMLE = @(par) -generateMLEfun(x, y, z, 1, 'STATZGATE', isoRef, par(1), par(2), par(3), par(4), par(5), par(6), par(7), vc1, vc2, vm);
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
        conRange95 = sqrt(chi2inv(0.95,Np)./diag(hessianMatrix));
        %                 chi2inv95 = chi2inv(0.95,Np);
        %                 conRange95 = sqrt(sqrt(diag(inv(chi2inv95*hessianMatrix)).^2));
        %         For DSL model
        %         Create empty sensitivity matrix
        %         sensitivityMatrix = zeros(length(x),Np);
        %         del = 1e-6;
        %         % Calculate sensitivity at every data point for each parameter
        %         for jj = 1:Np
        %             for kk = 1:length(x)
        %                 if jj == 1
        %                     sensitivityMatrix(kk,jj) = (computeStatZGATELoading(x(kk),y(kk),b01,delU1,beta,(1+del)*omega,kgate,cgate,gamma,vc1,vc2) - fitVals(kk))/(del*omega);
        %                 elseif jj == 2
        %                     sensitivityMatrix(kk,jj) = (computeStatZGATELoading(x(kk),y(kk),b01,delU1,(1+del)*beta,omega,kgate,cgate,gamma,vc1,vc2) - fitVals(kk))/(del*beta);
        %                 elseif jj == 3
        %                     sensitivityMatrix(kk,jj) = (computeStatZGATELoading(x(kk),y(kk),(1+del)*b01,delU1,beta,omega,kgate,cgate,gamma,vc1,vc2) - fitVals(kk))/(del*b01);
        %                 elseif jj == 4
        %                     sensitivityMatrix(kk,jj) = (computeStatZGATELoading(x(kk),y(kk),b01,(1+del)*delU1,beta,omega,kgate,cgate,gamma,vc1,vc2) - fitVals(kk))/(del*delU1);
        %                 elseif jj == 5
        %                     sensitivityMatrix(kk,jj) = (computeStatZGATELoading(x(kk),y(kk),b01,delU1,beta,omega,kgate,(1+del)*cgate,gamma,vc1,vc2) - fitVals(kk))/(del*cgate);
        %                 else
        %                     sensitivityMatrix(kk,jj) = (computeStatZGATELoading(x(kk),y(kk),b01,delU1,beta,omega,kgate,cgate,(1+del)*gamma,vc1,vc2) - fitVals(kk))/(del*gamma);
        %                 end
        %             end
        %         end
        %         % estimated Hessian Matrix for the data set (Non-linear parameter estimation
        %         % by Yonathan Bard (1974) pg. 178
        %         hessianMatrix = 1/stDevData^2*transpose(sensitivityMatrix)*sensitivityMatrix;
        %         % Confidence range given by chi squared distribution at Np degrees
        %         % of freedom (independent parameter conf intervals)
        %         conRange95 = sqrt(chi2inv(0.95,Np)./diag(hessianMatrix));
    case 'DSL'
        % Number of parameters
        Np = length(cell2mat(varargin(cell2mat(varargin)~=0)));
        % Calculate standard deviation of the data (not needed)
        stDevData = sqrt(1/(length(x)-Np) * sum((z-fitVals).^2));
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
                %                 chi2inv95 = chi2inv(0.95,length(cell2mat(varargin(cell2mat(varargin)~=0))));
                %                 if length(cell2mat(varargin(cell2mat(varargin)~=0))) == 6
                %                 else
                %                     hessianMatrix = [hessianMatrix(1,1) hessianMatrix(1,3) hessianMatrix(1,5);
                %                                      hessianMatrix(3,1) hessianMatrix(3,3) hessianMatrix(3,5);
                %                                      hessianMatrix(5,1) hessianMatrix(5,3) hessianMatrix(5,5)];
                %                 end
                %                     conRange95 = sqrt(sqrt(diag(inv(chi2inv95*hessianMatrix)).^2));
        end
    case 'DSL2'
        % Number of parameters
        Np = length(cell2mat(varargin(cell2mat(varargin)~=0)));
        % Calculate standard deviation of the data (not needed)
        stDevData = sqrt(1/(length(x)-Np) * sum((z-fitVals).^2));
        % Generate and solve global optimisation problem for confidence regions
        % based on isotherm model
        qs1a = varargin{1};
        qs2a = varargin{2};
        qs1b = varargin{3};
        qs2b = varargin{4};
        b01 = varargin{5};
        b02 = varargin{6};
        delU1 = varargin{7};
        delU2 = varargin{8};

        switch fittingMethod
            case 'MLE'
                % Create empty sensitivity matrix
                sensitivityMatrix = zeros(length(x),8);
                % Calculate sensitivity at every data point for each parameter
                for jj = 1:8
                    for kk = 1:length(x)
                        if jj == 1
                            sensitivityMatrix(kk,jj) = (((((1+del)*qs1a+qs1b./y(kk)))*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + (qs2a+qs2b./y(kk))*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))))) - fitVals(kk))/(del*qs1a);
                        elseif jj == 2
                            sensitivityMatrix(kk,jj) = (((qs1a+qs1b./y(kk))*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + ((1+del)*qs2a+qs2b./y(kk))*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))))) - fitVals(kk))/(del*qs2a);
                        elseif jj == 3
                            sensitivityMatrix(kk,jj) = ((((qs1a+(1+del)*qs1b./y(kk)))*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + (qs2a+qs2b./y(kk))*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))))) - fitVals(kk))/(del*qs1b);
                        elseif jj == 4
                            sensitivityMatrix(kk,jj) = (((qs1a+qs1b./y(kk))*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + (qs2a+(1+del).*qs2b./y(kk))*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))))) - fitVals(kk))/(del*qs2b);
                        elseif jj == 5
                            sensitivityMatrix(kk,jj) = (((qs1a+qs1b./y(kk))*(((1+del)*b01)*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(((1+del)*b01)*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + (qs2a+qs2b./y(kk))*(b02*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(b02*x(kk)*exp(delU2/(8.314*y(kk)))))) - fitVals(kk))/(del*b01);
                        elseif jj == 6
                            sensitivityMatrix(kk,jj) = (((qs1a+qs1b./y(kk))*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + (qs2a+qs2b./y(kk))*(((1+del)*b02)*x(kk)*exp(delU2/(8.314*y(kk))))/(1+(((1+del)*b02)*x(kk)*exp(delU2/(8.314*y(kk)))))) - fitVals(kk))/(del*b02);
                        elseif jj == 7
                            sensitivityMatrix(kk,jj) = (((qs1a+qs1b./y(kk))*(b01*x(kk)*exp((1+del)*delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp((1+del)*delU1/(8.314*y(kk))))) ...
                                + (qs2a+qs2b./y(kk))*((b02)*x(kk)*exp(delU2/(8.314*y(kk))))/(1+((b02)*x(kk)*exp(delU2/(8.314*y(kk)))))) - fitVals(kk))/(del*delU1);
                        else
                            sensitivityMatrix(kk,jj) = (((qs1a+qs1b./y(kk))*(b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))) ...
                                + (qs2a+qs2b./y(kk))*((b02)*x(kk)*exp((1+del)*delU2/(8.314*y(kk))))/(1+((b02)*x(kk)*exp((1+del)*delU2/(8.314*y(kk)))))) - fitVals(kk))/(del*delU2);
                        end
                    end
                end
                % estimated Hessian Matrix for the data set (Non-linear parameter estimation
                % by Yonathan Bard (1974) pg. 178
                hessianMatrix = 1/stDevData^2*transpose(sensitivityMatrix)*sensitivityMatrix;
                % Confidence range given by chi squared distribution at Np degrees
                % of freedom (independent parameter conf intervals)
                conRange95 = sqrt(chi2inv(0.95,Np)./diag(hessianMatrix));
                %                 chi2inv95 = chi2inv(0.95,Np);
                %                 conRange95 = sqrt(sqrt(diag(inv(chi2inv95*hessianMatrix)).^2));
            case 'WSS'
                Np = 8;
                Nt = length(x);
                dlogMLE = [];
                d2logMLE = [];
                deltaplus1mat = eye(Np).*(del);
                deltamat = eye(Np).*(del);
                partemp = [qs1a./isoRef(1), qs2a./isoRef(2),qs1b./isoRef(3), qs2b./isoRef(4),b01./isoRef(5),b02./isoRef(6),delU1./isoRef(7),delU2./isoRef(8)];
                logMLE = @(par) -generateMLEfun(x, y, z, 1, 'DSL2', isoRef, par(1), par(2), par(3), ...
                    par(4), par(5), par(6), par(7), par(8));

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
        stDevData = sqrt(1/(length(x)-Np) * sum((z-fitVals).^2));
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
                %                 chi2inv95 = chi2inv(0.95,Np);
                %                 conRange95 = sqrt(sqrt(diag(inv(chi2inv95*hessianMatrix)).^2));
        end
    case 'TSL'
        % Number of parameters
        Np = length(cell2mat(varargin(cell2mat(varargin)~=0)));
        % Calculate standard deviation of the data (not needed)
        stDevData = sqrt(1/(length(x)-Np) * sum((z-fitVals).^2));
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
                %                 chi2inv95 = chi2inv(0.95,Np);
                %                 conRange95 = sqrt(sqrt(diag(inv(chi2inv95*hessianMatrix)).^2));
        end
        % for DSS model
    case 'DSS'
        % Number of parameters
        Np = length(cell2mat(varargin(cell2mat(varargin)~=0)));
        % Calculate standard deviation of the data (not needed)
        stDevData = sqrt(1/(length(x)-Np) * sum((z-fitVals).^2));
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
                %                 chi2inv95 = chi2inv(0.95,length(cell2mat(varargin(cell2mat(varargin)~=0))));
                %                 if length(cell2mat(varargin(cell2mat(varargin)~=0))) == 7
                %                 else
                %                     hessianMatrix = [hessianMatrix(1,1) hessianMatrix(1,3) hessianMatrix(1,5) hessianMatrix(1,7);
                %                                      hessianMatrix(3,1) hessianMatrix(3,3) hessianMatrix(3,5) hessianMatrix(3,7);
                %                                      hessianMatrix(5,1) hessianMatrix(5,3) hessianMatrix(5,5) hessianMatrix(5,7);
                %                                      hessianMatrix(7,1) hessianMatrix(7,3) hessianMatrix(7,5) hessianMatrix(7,7)];
                %                 end
                %                 conRange95 = sqrt(sqrt(diag(inv(chi2inv95*hessianMatrix)).^2));
        end
    case 'DSS2'
        % Number of parameters
        Np = length(cell2mat(varargin(cell2mat(varargin)~=0)));
        % Calculate standard deviation of the data (not needed)
        stDevData = sqrt(1/(length(x)-Np) * sum((z-fitVals).^2));
        % Generate and solve global optimisation problem for confidence regions
        % based on isotherm model
        qs1 = varargin{1};
        qs2 = varargin{2};
        b01 = varargin{3};
        b02 = varargin{4};
        delU1 = varargin{5};
        delU2 = varargin{6};
        gamma = varargin{7};

        Np = 7;
        Nt = length(x);
        dlogMLE = [];
        d2logMLE = [];
        deltaplus1mat = eye(Np).*(del);
        deltamat = eye(Np).*(del);
        partemp = [qs1./isoRef(1), qs2./isoRef(2),b01./isoRef(3),b02./isoRef(4),delU1./isoRef(5),delU2./isoRef(6),...
            gamma./isoRef(7)];
        logMLE = @(par) -generateMLEfun(x, y, z, 1, 'DSS2', isoRef, par(1), par(2), par(3), ...
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
        %                 chi2inv95 = chi2inv(0.95,length(cell2mat(varargin(cell2mat(varargin)~=0))));
        %                 if length(cell2mat(varargin(cell2mat(varargin)~=0))) == 7
        %                 else
        %                     hessianMatrix = [hessianMatrix(1,1) hessianMatrix(1,3) hessianMatrix(1,5) hessianMatrix(1,7);
        %                                      hessianMatrix(3,1) hessianMatrix(3,3) hessianMatrix(3,5) hessianMatrix(3,7);
        %                                      hessianMatrix(5,1) hessianMatrix(5,3) hessianMatrix(5,5) hessianMatrix(5,7);
        %                                      hessianMatrix(7,1) hessianMatrix(7,3) hessianMatrix(7,5) hessianMatrix(7,7)];
        %                 end
        %                 conRange95 = sqrt(sqrt(diag(inv(chi2inv95*hessianMatrix)).^2));

        % for DSS model
    case 'TOTH'
        % Number of parameters
        Np = length(cell2mat(varargin(cell2mat(varargin)~=0)));
        % Calculate standard deviation of the data (not needed)
        stDevData = sqrt(1/(length(x)-Np) * sum((z-fitVals).^2));
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
            case 'WSS'
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
                %                 chi2inv95 = chi2inv(0.95,length(cell2mat(varargin(cell2mat(varargin)~=0))));
                %                 if length(cell2mat(varargin(cell2mat(varargin)~=0))) == 7
                %                 else
                %                     hessianMatrix = [hessianMatrix(1,1) hessianMatrix(1,3) hessianMatrix(1,5) hessianMatrix(1,7);
                %                                      hessianMatrix(3,1) hessianMatrix(3,3) hessianMatrix(3,5) hessianMatrix(3,7);
                %                                      hessianMatrix(5,1) hessianMatrix(5,3) hessianMatrix(5,5) hessianMatrix(5,7);
                %                                      hessianMatrix(7,1) hessianMatrix(7,3) hessianMatrix(7,5) hessianMatrix(7,7)];
                %                 end
                %                 conRange95 = sqrt(sqrt(diag(inv(chi2inv95*hessianMatrix)).^2));
            case 'MLE'
                Np = 7;
                Nt = length(x);
                dlogMLE = [];
                d2logMLE = [];
                deltaplus1mat = eye(Np).*(del);
                deltamat = eye(Np).*(del);
                partemp = [qs1./isoRef(1), 0, b01./isoRef(3), 0, delU1./isoRef(5), 0, toth./isoRef(7)];
                logMLE = @(par) -generateMLEfun(x, y, z, 1, 'TOTH', isoRef, par(1), par(2), par(3), ...
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
                conRange95 = sqrt(chi2inv(0.95, length(cell2mat(varargin(cell2mat(varargin)~=0))))./diag(hessianMatrix));
                %                 chi2inv95 = chi2inv(0.95,length(cell2mat(varargin(cell2mat(varargin)~=0))));
                %                 if length(cell2mat(varargin(cell2mat(varargin)~=0))) == 7
                %                 else
                %                     hessianMatrix = [hessianMatrix(1,1) hessianMatrix(1,3) hessianMatrix(1,5) hessianMatrix(1,7);
                %                                      hessianMatrix(3,1) hessianMatrix(3,3) hessianMatrix(3,5) hessianMatrix(3,7);
                %                                      hessianMatrix(5,1) hessianMatrix(5,3) hessianMatrix(5,5) hessianMatrix(5,7);
                %                                      hessianMatrix(7,1) hessianMatrix(7,3) hessianMatrix(7,5) hessianMatrix(7,7)];
                %                 end
                %                 conRange95 = sqrt(sqrt(diag(inv(chi2inv95*hessianMatrix)).^2));
        end
    case 'TOTH2'
        % Number of parameters
        Np = length(cell2mat(varargin(cell2mat(varargin)~=0)));
        % Calculate standard deviation of the data (not needed)
        stDevData = sqrt(1/(length(x)-Np) * sum((z-fitVals).^2));
        % Generate and solve global optimisation problem for confidence regions
        % based on isotherm model
        qs1 = varargin{1};
        qs2 = varargin{2};
        b01 = varargin{3};
        b02 = varargin{4};
        delU1 = varargin{5};
        delU2 = varargin{6};
        toth0 = varargin{7};
        totha = varargin{8};
        switch fittingMethod
            case 'MLE'
                % Create empty sensitivity matrix
                sensitivityMatrix = zeros(length(x),8);
                % Calculate sensitivity at every data point for each parameter
                for jj = 1:8
                    for kk = 1:length(x)
                        if jj == 1
                            sensitivityMatrix(kk,jj) = abs(((((1+del)*qs1)*b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^(toth0 + totha.*(1-298.15/y(kk)))).^(1./(toth0 + totha.*(1-298.15/y(kk)))) ...
                                - fitVals(kk))/(del*qs1));
                        elseif jj == 2
                            sensitivityMatrix(kk,jj) = abs((((1+del)*qs1*b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^(toth0 + totha.*(1-298.15/y(kk)))).^(1./(toth0 + totha.*(1-298.15/y(kk)))) ...
                                - fitVals(kk))/(del*qs1));
                        elseif jj == 3
                            sensitivityMatrix(kk,jj) = abs(((qs1*(1+del)*b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^(toth0 + totha.*(1-298.15/y(kk)))).^(1./(toth0 + totha.*(1-298.15/y(kk)))) ...
                                - fitVals(kk))/(del*b01));
                        elseif jj == 4
                            sensitivityMatrix(kk,jj) = abs(((qs1*(1+del)*b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^(toth0 + totha.*(1-298.15/y(kk)))).^(1./(toth0 + totha.*(1-298.15/y(kk)))) ...
                                - fitVals(kk))/(del*b01));
                        elseif jj == 5
                            sensitivityMatrix(kk,jj) = abs(((qs1*b01*x(kk)*exp((1+del)*delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^(toth0 + totha.*(1-298.15/y(kk)))).^(1./(toth0 + totha.*(1-298.15/y(kk)))) ...
                                - fitVals(kk))/(del*delU1));
                        elseif jj == 6
                            sensitivityMatrix(kk,jj) = abs(((qs1*b01*x(kk)*exp((1+del)*delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^(toth0 + totha.*(1-298.15/y(kk)))).^(1./(toth0 + totha.*(1-298.15/y(kk)))) ...
                                - fitVals(kk))/(del*delU1));
                        elseif jj == 7
                            sensitivityMatrix(kk,jj) = abs(((qs1*b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^((1+del)*toth0 + totha.*(1-298.15/y(kk)))).^(1./((1+del)*toth0 + totha.*(1-298.15/y(kk)))) ...
                                - fitVals(kk))/(del*toth0));
                        else
                            sensitivityMatrix(kk,jj) = abs(((qs1*b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^(toth0 + (1+del)*totha.*(1-298.15/y(kk)))).^(1./(toth0 + (1+del)*totha.*(1-298.15/y(kk)))) ...
                                - fitVals(kk))/(del*totha));
                        end
                    end
                end
                % Hessian Matrix for the data set (Non-linear parameter estimation
                % by Yonathan Bard (1974) pg. 178)
                hessianMatrix = 1/stDevData^2*transpose(sensitivityMatrix)*sensitivityMatrix;
                % Confidence range given by chi squared distribution at Np degrees
                % of freedom (independent parameter conf intervals)
                conRange95 = sqrt(chi2inv(0.95,Np)./diag(hessianMatrix));
                %                 hessianMatrix = [hessianMatrix(1,1) hessianMatrix(1,3) hessianMatrix(1,5) hessianMatrix(1,7) hessianMatrix(1,8);
                %                     hessianMatrix(3,1) hessianMatrix(3,3) hessianMatrix(3,5) hessianMatrix(3,7) hessianMatrix(3,8);
                %                     hessianMatrix(5,1) hessianMatrix(5,3) hessianMatrix(5,5) hessianMatrix(5,7) hessianMatrix(5,8);
                %                     hessianMatrix(7,1) hessianMatrix(7,3) hessianMatrix(7,5) hessianMatrix(7,7) hessianMatrix(7,8);
                %                     hessianMatrix(8,1) hessianMatrix(8,3) hessianMatrix(8,5) hessianMatrix(8,7) hessianMatrix(8,8)];
                %                 chi2inv95 = chi2inv(0.95,Np);
                %                 conRange95 = sqrt(sqrt(diag(inv(chi2inv95*hessianMatrix)).^2));
            case 'WSS'
                Np = 8;
                Nt = length(x);
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
                %                 chi2inv95 = chi2inv(0.95,Np);
                %                 conRange95 = sqrt(sqrt(diag(inv(chi2inv95*hessianMatrix)).^2));
        end
    case 'TOTH3'
        % Number of parameters
        Np = 6;
        % Calculate standard deviation of the data (not needed)
        stDevData = sqrt(1/(length(x)-Np) * sum((z-fitVals).^2));
        % Generate and solve global optimisation problem for confidence regions
        % based on isotherm model
        qs10 = varargin{1};
        qs2 = varargin{2};
        b01 = varargin{3};
        b02 = varargin{4};
        delU1 = varargin{5};
        delU2 = varargin{6};
        toth0 = varargin{7};
        totha = varargin{8};
        chi = varargin{9};
        switch fittingMethod
            case 'MLE'
                % Create empty sensitivity matrix
                sensitivityMatrix = zeros(length(x),6);
                % Calculate sensitivity at every data point for each parameter
                for jj = 1:6
                    for kk = 1:length(x)
                        if jj == 1
                            sensitivityMatrix(kk,jj) = ((((1+del)*qs10.*exp(chi*(1-y(kk)./298.15)))*b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^(toth0 + totha.*(1-298.15/y(kk)))).^(1./(toth0 + totha.*(1-298.15/y(kk)))) ...
                                - fitVals(kk))/(del*qs10);
                        elseif jj == 2
                            sensitivityMatrix(kk,jj) = ((qs10.*exp(chi*(1-y(kk)./298.15))*(1+del)*b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+((1+del)*b01*x(kk)*exp(delU1/(8.314*y(kk))))^(toth0 + totha.*(1-298.15/y(kk)))).^(1./(toth0 + totha.*(1-298.15/y(kk)))) ...
                                - fitVals(kk))/(del*b01);
                        elseif jj == 3
                            sensitivityMatrix(kk,jj) = ((qs10.*exp(chi*(1-y(kk)./298.15))*b01*x(kk)*exp((1+del)*delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp((1+del)*delU1/(8.314*y(kk))))^(toth0 + totha.*(1-298.15/y(kk)))).^(1./(toth0 + totha.*(1-298.15/y(kk)))) ...
                                - fitVals(kk))/(del*delU1);
                        elseif jj == 4
                            sensitivityMatrix(kk,jj) = ((qs10.*exp(chi*(1-y(kk)./298.15))*b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^((1+del)*toth0 + totha.*(1-298.15/y(kk)))).^(1./((1+del)*toth0 + totha.*(1-298.15/y(kk)))) ...
                                - fitVals(kk))/(del*toth0);
                        elseif jj == 5
                            sensitivityMatrix(kk,jj) = ((qs10.*exp(chi*(1-y(kk)./298.15))*b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^(toth0 + (1+del)*totha.*(1-298.15/y(kk)))).^(1./(toth0 + (1+del)*totha.*(1-298.15/y(kk)))) ...
                                - fitVals(kk))/(del*totha);
                        else
                            sensitivityMatrix(kk,jj) = ((qs10.*exp((1+del)*chi*(1-y(kk)./298.15))*b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^(toth0 + totha.*(1-298.15/y(kk)))).^(1./(toth0 + totha.*(1-298.15/y(kk)))) ...
                                - fitVals(kk))/(del*chi);
                        end
                    end
                end
                % Hessian Matrix for the data set (Non-linear parameter estimation
                % by Yonathan Bard (1974) pg. 178)
                hessianMatrix = 1/stDevData^2*transpose(sensitivityMatrix)*sensitivityMatrix;
                % Confidence range given by chi squared distribution at Np degrees
                % of freedom (independent parameter conf intervals)
                conRange95 = sqrt(chi2inv(0.95,Np)./diag(hessianMatrix));
                %                 chi2inv95 = chi2inv(0.95,Np);
                %                 conRange95 = sqrt(sqrt(diag(inv(chi2inv95*hessianMatrix)).^2));
            case 'WSS'
                Np = 7;
                Nt = length(x);
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
                chi2inv95 = chi2inv(0.95,Np);
                conRange95 = sqrt(sqrt(diag(inv(chi2inv95*hessianMatrix)).^2));
        end
    case 'TOTHCHEM'
        % Number of parameters
        Np = 10;
        % Calculate standard deviation of the data (not needed)
        stDevData = sqrt(1/(length(x)-Np) * sum((z-fitVals).^2));
        % Generate and solve global optimisation problem for confidence regions
        % based on isotherm model
        qs10 = varargin{1};
        qs2 = varargin{2};
        b01 = varargin{3};
        b02 = varargin{4};
        delU1 = varargin{5};
        delU2 = varargin{6};
        toth0 = varargin{7};
        totha = varargin{8};
        chi = varargin{9};
        qsC = varargin{10};
        b0C = varargin{11};
        delUC = varargin{12};
        EaC = varargin{13};
        switch fittingMethod
            case 'MLE'
                % Create empty sensitivity matrix
                sensitivityMatrix = zeros(length(x),10);
                % Calculate sensitivity at every data point for each parameter
                for jj = 1:Np
                    for kk = 1:length(x)
                        if jj == 1
                            sensitivityMatrix(kk,jj) = ((((1+del)*qs10.*exp(chi*(1-y(kk)./298.15)))*b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^(toth0 + totha.*(1-298.15/y(kk)))).^(1./(toth0 + totha.*(1-298.15/y(kk)))) ...
                                + exp(-EaC/(8.314.*y(kk)))*qsC.*b0C.*x(kk).*exp(delUC./(8.314.*y(kk)))./(1+b0C.*x(kk).*exp(delUC./(8.314.*y(kk)))) ...
                                - fitVals(kk))/(del*qs10);
                        elseif jj == 2
                            sensitivityMatrix(kk,jj) = ((qs10.*exp(chi*(1-y(kk)./298.15))*(1+del)*b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+((1+del)*b01*x(kk)*exp(delU1/(8.314*y(kk))))^(toth0 + totha.*(1-298.15/y(kk)))).^(1./(toth0 + totha.*(1-298.15/y(kk)))) ...
                                + exp(-EaC/(8.314.*y(kk)))*qsC.*b0C.*x(kk).*exp(delUC./(8.314.*y(kk)))./(1+b0C.*x(kk).*exp(delUC./(8.314.*y(kk)))) ...
                                - fitVals(kk))/(del*b01);
                        elseif jj == 3
                            sensitivityMatrix(kk,jj) = ((qs10.*exp(chi*(1-y(kk)./298.15))*b01*x(kk)*exp((1+del)*delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp((1+del)*delU1/(8.314*y(kk))))^(toth0 + totha.*(1-298.15/y(kk)))).^(1./(toth0 + totha.*(1-298.15/y(kk)))) ...
                                + exp(-EaC/(8.314.*y(kk)))*qsC.*b0C.*x(kk).*exp(delUC./(8.314.*y(kk)))./(1+b0C.*x(kk).*exp(delUC./(8.314.*y(kk)))) ...
                                - fitVals(kk))/(del*delU1);
                        elseif jj == 4
                            sensitivityMatrix(kk,jj) = ((qs10.*exp(chi*(1-y(kk)./298.15))*b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^((1+del)*toth0 + totha.*(1-298.15/y(kk)))).^(1./((1+del)*toth0 + totha.*(1-298.15/y(kk)))) ...
                                + exp(-EaC/(8.314.*y(kk)))*qsC.*b0C.*x(kk).*exp(delUC./(8.314.*y(kk)))./(1+b0C.*x(kk).*exp(delUC./(8.314.*y(kk)))) ...
                                - fitVals(kk))/(del*toth0);
                        elseif jj == 5
                            sensitivityMatrix(kk,jj) = ((qs10.*exp(chi*(1-y(kk)./298.15))*b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^(toth0 + (1+del)*totha.*(1-298.15/y(kk)))).^(1./(toth0 + (1+del)*totha.*(1-298.15/y(kk)))) ...
                                + exp(-EaC/(8.314.*y(kk)))*qsC.*b0C.*x(kk).*exp(delUC./(8.314.*y(kk)))./(1+b0C.*x(kk).*exp(delUC./(8.314.*y(kk)))) ...
                                - fitVals(kk))/(del*totha);
                        elseif jj == 6
                            sensitivityMatrix(kk,jj) = ((qs10.*exp((1+del)*chi*(1-y(kk)./298.15))*b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^(toth0 + totha.*(1-298.15/y(kk)))).^(1./(toth0 + totha.*(1-298.15/y(kk)))) ...
                                + exp(-EaC/(8.314.*y(kk)))*qsC.*b0C.*x(kk).*exp(delUC./(8.314.*y(kk)))./(1+b0C.*x(kk).*exp(delUC./(8.314.*y(kk)))) ...
                                - fitVals(kk))/(del*chi);
                        elseif jj == 7
                            sensitivityMatrix(kk,jj) = ((qs10.*exp(chi*(1-y(kk)./298.15))*b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^(toth0 + totha.*(1-298.15/y(kk)))).^(1./(toth0 + totha.*(1-298.15/y(kk)))) ...
                                + exp(-EaC/(8.314.*y(kk)))*(1+del)*qsC.*b0C.*x(kk).*exp(delUC./(8.314.*y(kk)))./(1+b0C.*x(kk).*exp(delUC./(8.314.*y(kk)))) ...
                                - fitVals(kk))/(del*qsC);
                        elseif jj == 8
                            sensitivityMatrix(kk,jj) = ((qs10.*exp(chi*(1-y(kk)./298.15))*b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^(toth0 + totha.*(1-298.15/y(kk)))).^(1./(toth0 + totha.*(1-298.15/y(kk)))) ...
                                + exp(-EaC/(8.314.*y(kk)))*qsC.*(1+del).*b0C.*x(kk).*exp(delUC./(8.314.*y(kk)))./(1+(1+del).*b0C.*x(kk).*exp(delUC./(8.314.*y(kk)))) ...
                                - fitVals(kk))/(del*b0C);
                        elseif jj == 9
                            sensitivityMatrix(kk,jj) = ((qs10.*exp(chi*(1-y(kk)./298.15))*b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^(toth0 + totha.*(1-298.15/y(kk)))).^(1./(toth0 + totha.*(1-298.15/y(kk)))) ...
                                + exp(-EaC/(8.314.*y(kk)))*qsC.*b0C.*x(kk).*exp((1+del)*delUC./(8.314.*y(kk)))./(1+b0C.*x(kk).*exp((1+del)*delUC./(8.314.*y(kk)))) ...
                                - fitVals(kk))/(del*delUC);
                        else
                            sensitivityMatrix(kk,jj) = ((qs10.*exp(chi*(1-y(kk)./298.15))*b01*x(kk)*exp(delU1/(8.314*y(kk))))/(1+(b01*x(kk)*exp(delU1/(8.314*y(kk))))^(toth0 + totha.*(1-298.15/y(kk)))).^(1./(toth0 + totha.*(1-298.15/y(kk)))) ...
                                + exp(-(1+del)*EaC/(8.314.*y(kk)))*qsC.*b0C.*x(kk).*exp(delUC./(8.314.*y(kk)))./(1+b0C.*x(kk).*exp(delUC./(8.314.*y(kk)))) ...
                                - fitVals(kk))/(del*EaC);
                        end
                    end
                end
                % Hessian Matrix for the data set (Non-linear parameter estimation
                % by Yonathan Bard (1974) pg. 178)
                hessianMatrix = 1/stDevData^2*transpose(sensitivityMatrix)*sensitivityMatrix;
                % Confidence range given by chi squared distribution at Np degrees
                % of freedom (independent parameter conf intervals)
                conRange95 = sqrt(chi2inv(0.95,Np)./diag(hessianMatrix));
                %                 chi2inv95 = chi2inv(0.95,Np);
                %                 conRange95 = sqrt(sqrt(diag(inv(chi2inv95*hessianMatrix)).^2))
        end
        % for Virial model
    case 'VIRIAL'
        % Number of parameters
        Np = length(cell2mat(varargin(cell2mat(varargin)~=0)));
        % Calculate standard deviation of the data (not needed)
        stDevData = sqrt(1/(length(x)-Np) * sum((log(x)-fitVals').^2));
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
            case 'MLE'
                % Create empty sensitivity matrix
                sensitivityMatrix = zeros(length(x),8);
                % Calculate sensitivity at every data point for each parameter
                for jj = 1:Np
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
                hessianMatrix = hessianMatrix(1:Np,1:Np);
                conRange95 = sqrt(chi2inv(0.95,Np)./diag(hessianMatrix));
                chi2inv95 = chi2inv(0.95,Np);
                conRange95 = sqrt(sqrt(diag(inv(chi2inv95*hessianMatrix)).^2));
            case 'WSS'
                Nt = length(x);
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
                %                 chi2inv95 = chi2inv(0.95,Np);
                %                 conRange95 = sqrt(sqrt(diag(inv(chi2inv95*hessianMatrix)).^2));
        end
    case 'VIRIAL2'
        % Number of parameters
        Np = length(cell2mat(varargin));
        % Calculate standard deviation of the data (not needed)
        stDevData = sqrt(1/(length(x)-Np) * sum((log(x)-fitVals').^2));
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
            case 'MLE'
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
                %                 chi2inv95 = chi2inv(0.95,Np);
                %                 conRange95 = sqrt(sqrt(diag(inv(chi2inv95*hessianMatrix)).^2));
            case 'WSS'
                Nt = length(x);
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
                %                 chi2inv95 = chi2inv(0.95,Np);
                %                 conRange95 = sqrt(sqrt(diag(inv(chi2inv95*hessianMatrix)).^2));
        end
end
end