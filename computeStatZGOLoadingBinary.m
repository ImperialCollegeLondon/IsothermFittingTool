function [qa, qb, qT] = computeStatZGOLoadingBinary(P,T,b01A,delU1A,delU2A,kgateA,cgateA,betaA,omegaA,b01B,delU1B,delU2B,kgateB,cgateB,betaB,omegaB,vc,yA)

pA = P.*yA;
pB = P.*(1-yA);

delUA = (delU1A+delU2A)./2 + (delU2A-delU1A)./2.*tanh(kgateA.*b01A.*exp(delU1A./(8.314.*T)).*pA - cgateA);
delUB = (delU1B+delU2B)./2 + (delU2B-delU1B)./2.*tanh(kgateB.*b01B.*exp(delU1B./(8.314.*T)).*pB - cgateB);

delUA = (delU1A+delU2A)./2 + (delU2A-delU1A)./2.*tanh((kgateA.*(b01A.*exp(delU1A./(8.314.*T)).*pA)./(1+b01A.*exp(delU1A./(8.314.*T)).*pA) - cgateA)./2);
delUB = (delU1B+delU2B)./2 + (delU2B-delU1B)./2.*tanh((kgateB.*(b01B.*exp(delU1B./(8.314.*T)).*pB)./(1+b01B.*exp(delU1B./(8.314.*T)).*pB) - cgateB)./2);
% delUB = (delU1B+delU2B)./2 + (delU2B-delU1B)./2.*tanh(kgateA.*(b01A.*exp(delU1A./(8.314.*T)).*pA)./(1+b01A.*exp(delU1A./(8.314.*T)).*pA) - cgateA);


bA = b01A.*exp(delUA./(8.314.*T));
bB = b01B.*exp(delUB./(8.314.*T));



isothermNumeratorA = 0;
isothermNumeratorB = 0;
isothermDenominator = 0;

% evaluate double summation of over integers for i (jj) and j (kk)
for kk = 0:1000
    for jj = 0:1000
        % subject to constraints
        % if (jj*betaA + kk*betaB) <= vc && kk+jj >= 2
        if (jj*betaA + kk*betaB) <= max(omegaA.*betaA,omegaB.*betaB) && kk+jj >= 2
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