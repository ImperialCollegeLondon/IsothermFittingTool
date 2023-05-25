function qa = computeStatZGATELoading(P,T,b01,delU1,beta,omega,kgate,cgate,gamma,vc1,vc2)

b = b01.*exp(delU1./(8.314.*T));
vc = (vc1+vc2)./2 + (vc2-vc1)./2.*tanh(kgate.*P - cgate);
isothermNumerator = 0;
isothermDenominator =0;

for jj = 1:omega
    isothermNumerator = isothermNumerator + ((b.*P).^jj./(factorial(jj-1))).* ...
        ((1-jj*beta./vc)./(1-beta./vc)).^jj;
    isothermDenominator = isothermDenominator + ((b.*P).^jj./(factorial(jj))).* ...
        ((1-jj*beta./vc)./(1-beta./vc)).^jj;

%         isothermNumerator = isothermNumerator + (b.*P).^(kgate+jj)./(factorial(jj-1)).* ...
%         ((1-jj*beta./vc)./(1-beta./vc)).^jj;
%     isothermDenominator = isothermDenominator + (b.*P).^(kgate+jj)./(factorial(jj)).* ...
%         ((1-jj*beta./vc)./(1-beta./vc)).^jj;
end
qa = isothermNumerator./(1+isothermDenominator.^gamma).^(1/gamma);
% qa = (isothermNumerator./(1+isothermDenominator)).^(1./gamma);

end