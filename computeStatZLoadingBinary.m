function [qa, qb, qT] = computeStatZLoadingBinary(P,T,b01A,delU1A,betaA,omegaA,b01B,delU1B,betaB,omegaB,vc,yA)
bA = b01A.*exp(delU1A./(8.314.*T));
bB = b01B.*exp(delU1B./(8.314.*T));

pA = P.*yA;
pB = P.*(1-yA);

isothermNumeratorA = 0;
isothermNumeratorB = 0;
isothermDenominator = 0;

% evaluate double summation of over integers for i (jj) and j (kk)
for kk = 0:1000
    for jj = 0:1000
        % subject to constraints
        if (jj*betaA + kk*betaB) <= vc && kk+jj >= 2
            if     kk == 0
                % Compute numnerator and denominator amount for A for
                % non-competitive case
                isothermNumeratorA =   isothermNumeratorA   + (bA.*pA).^jj./factorial(jj-1).*(1-jj*betaA./vc).^jj;
                isothermDenominator =  isothermDenominator  + (bA.*pA).^jj./factorial(jj).*  (bB.*pB).^kk./factorial(kk).*(1-jj.*betaA./vc-kk.*betaB./vc).^(jj+kk);
            elseif jj == 0
                % Compute numnerator and denominator amount for B for
                % non-competitive case
                isothermNumeratorB =   isothermNumeratorB   + (bB.*pB).^kk./factorial(kk-1).*(1-kk*betaB./vc).^kk;
                isothermDenominator =  isothermDenominator  + (bB.*pB).^kk./factorial(kk).*  (bA.*pA).^jj./factorial(jj).*(1-jj.*betaA./vc-kk.*betaB./vc).^(jj+kk);
            else
                % Compute numnerator and denominator amount for A for
                % competitive case
                isothermNumeratorA =   isothermNumeratorA   + (bA.*pA).^jj./factorial(jj-1).*(bB.*pB).^kk./factorial(kk).*(1-jj.*betaA./vc-kk.*betaB./vc).^(jj+kk);
                % Compute numnerator and denominator amount for B for
                % competitive  case
                isothermNumeratorB =   isothermNumeratorB   + (bB.*pB).^kk./factorial(kk-1).*(bA.*pA).^jj./factorial(jj).*(1-jj.*betaA./vc-kk.*betaB./vc).^(jj+kk);
                isothermDenominator =  isothermDenominator  + (bB.*pB).^kk./factorial(kk).*  (bA.*pA).^jj./factorial(jj).*(1-jj.*betaA./vc-kk.*betaB./vc).^(jj+kk);
            end
        end
    end
end

qa = ((bA.*pA+isothermNumeratorA)./(1+bA.*pA+bB.*pB+isothermDenominator))';
qb = ((bB.*pB+isothermNumeratorB)./(1+bA.*pA+bB.*pB+isothermDenominator))';
qT = qa+qb;
end