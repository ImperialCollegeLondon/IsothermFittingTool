function [qa, qb, qT] = computeStatSTALoadingBinary(P,T,b01A,b02A,delU1A,delU2A,betaA,b01B,b02B,delU1B,delU2B,betaB,kgate,cgate,delvc,omegaA,omegaB,vc,yA,varargin)

if length(varargin) ~= 0
    gammaA = varargin{1};
    gammaB = varargin{2};
else
    gammaA = 1;
    gammaB = 1;
end

% delvc2 = delvc.*vc;

betavdWA = vc./omegaA;
betavdWB = vc./omegaB;

% betaA = betavdWA;
% betaB = betavdWB;


vc1 = vc - delvc.*vc;

omega1A = floor(vc1./betavdWA);
omega1B = floor(vc1./betavdWB);



pA = P.*yA;
pB = P.*(1-yA);


b1A = b01A.*exp(delU1A./(8.314.*T));
b2A = b02A.*exp(delU2A./(8.314.*T));

b1B = b01B.*exp(delU1B./(8.314.*T));
b2B = b02B.*exp(delU2B./(8.314.*T));

isothermDenominator1 = 0;
isothermDenominator2 = 0;
isothermNumerator1A   = 0;
isothermNumerator2A   = 0;
isothermNumerator1B   = 0;
isothermNumerator2B   = 0;


if omega1A>0
    % evaluate double summation of over integers for i (jj) and j (kk) FOR
    % PHASE 1
    for kk = 0:1000
        for jj = 0:1000
            % subject to constraints
            if (jj*betaA + kk*betaB) <= min(omega1A.*betavdWA,omega1B.*betavdWB) && kk+jj >= 2
                if     kk == 0
                    % Compute numnerator and denominator amount for A for
                    % non-competitive case
                    isothermNumerator1A =   isothermNumerator1A   + ((b1A.*pA).^gammaA).^jj./factorial(jj-1).*(1-jj*betaA./vc1).^jj;
                    isothermDenominator1 =  isothermDenominator1  + ((b1A.*pA).^gammaA).^jj./factorial(jj).*  ((b1B.*pB).^gammaB).^kk./factorial(kk).*(1-jj.*betaA./vc1-kk.*betaB./vc1).^(jj+kk);
                elseif jj == 0
                    % Compute numnerator and denominator amount for B for
                    % non-competitive case
                    isothermNumerator1B =   isothermNumerator1B   + ((b1B.*pB).^gammaB).^kk./factorial(kk-1).*(1-kk*betaB./vc1).^kk;
                    isothermDenominator1 =  isothermDenominator1  + ((b1B.*pB).^gammaB).^kk./factorial(kk).*  ((b1A.*pA).^gammaA).^jj./factorial(jj).*(1-jj.*betaA./vc1-kk.*betaB./vc1).^(jj+kk);
                else
                    % Compute numnerator and denominator amount for A for
                    % competitive case
                    isothermNumerator1A =   isothermNumerator1A   + ((b1A.*pA).^gammaA).^jj./factorial(jj-1).*((b1B.*pB).^gammaB).^kk./factorial(kk).*(1-jj.*betaA./vc1-kk.*betaB./vc1).^(jj+kk);
                    % Compute numnerator and denominator amount for B for
                    % competitive  case
                    isothermNumerator1B =   isothermNumerator1B   + ((b1B.*pB).^gammaB).^kk./factorial(kk-1).*((b1A.*pA).^gammaA).^jj./factorial(jj).*(1-jj.*betaA./vc1-kk.*betaB./vc1).^(jj+kk);
                    isothermDenominator1 =  isothermDenominator1  + ((b1B.*pB).^gammaB).^kk./factorial(kk).*  ((b1A.*pA).^gammaA).^jj./factorial(jj).*(1-jj.*betaA./vc1-kk.*betaB./vc1).^(jj+kk);
                end
            end
        end
    end
end

% evaluate double summation of over integers for i (jj) and j (kk) FOR
% PHASE 2
for kk = 0:1000
    for jj = 0:1000
        % subject to constraints
        if (jj*betaA + kk*betaB) <= min(omegaA.*betavdWA,omegaB.*betavdWB) && kk+jj >= 2
            if  kk == 0
                % Compute numnerator and denominator amount for A for
                % non-competitive case
                isothermNumerator2A =   isothermNumerator2A   + ((b2A.*pA).^gammaA).^jj./factorial(jj-1).*(1-jj*betaA./vc).^jj;
                isothermDenominator2 =  isothermDenominator2  + ((b2A.*pA).^gammaA).^jj./factorial(jj).*  ((b2B.*pB).^gammaB).^kk./factorial(kk).*(1-jj.*betaA./vc-kk.*betaB./vc).^(jj+kk);
            elseif jj == 0
                % Compute numnerator and denominator amount for B for
                % non-competitive case
                isothermNumerator2B =   isothermNumerator2B   + ((b2B.*pB).^gammaB).^kk./factorial(kk-1).*(1-kk*betaB./vc1).^kk;
                isothermDenominator2 =  isothermDenominator2  + ((b2B.*pB).^gammaB).^kk./factorial(kk).*  ((b2A.*pA).^gammaA).^jj./factorial(jj).*(1-jj.*betaA./vc-kk.*betaB./vc).^(jj+kk);
            else
                % Compute numnerator and denominator amount for A for
                % competitive case
                isothermNumerator2A =   isothermNumerator2A   + ((b2A.*pA).^gammaA).^jj./factorial(jj-1).*((b2B.*pB).^gammaB).^kk./factorial(kk).*(1-jj.*betaA./vc-kk.*betaB./vc).^(jj+kk);
                % Compute numnerator and denominator amount for B for
                % competitive  case
                isothermNumerator2B =   isothermNumerator2B   + ((b2B.*pB).^gammaB).^kk./factorial(kk-1).*((b2A.*pA).^gammaA).^jj./factorial(jj).*(1-jj.*betaA./vc-kk.*betaB./vc).^(jj+kk);
                isothermDenominator2 =  isothermDenominator2  + ((b2B.*pB).^gammaB).^kk./factorial(kk).*  ((b2A.*pA).^gammaA).^jj./factorial(jj).*(1-jj.*betaA./vc-kk.*betaB./vc).^(jj+kk);
            end
        end
    end
end


qa1 = (((b1A.*pA).^gammaA+isothermNumerator1A)./(1+(b1A.*pA).^gammaA+(b1B.*pB).^gammaB+isothermDenominator1))';
qa2 = (((b2A.*pA).^gammaA+isothermNumerator2A)./(1+(b2A.*pA).^gammaA+(b2B.*pB).^gammaB+isothermDenominator2))';

qb1 = (((b1B.*pB).^gammaB+isothermNumerator1B)./(1+(b1A.*pA).^gammaA+(b1B.*pB).^gammaB+isothermDenominator1))';
qb2 = (((b2B.*pB).^gammaB+isothermNumerator2B)./(1+(b2A.*pA).^gammaA+(b2B.*pB).^gammaB+isothermDenominator2))';

 
yval = ((1+(b2A.*pA).^gammaA+(b2B.*pB).^gammaB+isothermDenominator2))./((1+(b1A.*pA).^gammaA+(b1B.*pB).^gammaB+isothermDenominator1)).*exp(-(kgate.*(8.314) - T.*cgate.*(8.314))./(T.*(8.314)));
 
sigmaval = yval./(1+yval);

qa = (1-sigmaval').*qa1 + sigmaval'.*qa2;
qb = (1-sigmaval').*qb1 + sigmaval'.*qb2;
qT = qa+qb;


% subplot(1,2,1)
% hold on
% plot(P, (isothermNumerator1))
% plot(P, (isothermNumerator2).*exp(-(kgate./T - cgate)))
% set(gca,'YScale','log')
% subplot(1,2,2)
% hold on
% plot(P, (((1+(b1A.*pA).^gammaA+(b1B.*pB).^gammaB+isothermDenominator1))))
% hold on
% plot(P, ((1+(b2A.*pA).^gammaA+(b2B.*pB).^gammaB+isothermDenominator2)).*exp(-(kgate./T - cgate)))
% set(gca,'YScale','log')

end