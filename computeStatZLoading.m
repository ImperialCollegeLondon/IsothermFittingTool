function qa = computeStatZLoading(P,T,b01,delU1,beta,omega,vc)

b = b01.*exp(delU1./(8.314.*T));

isothermNumerator = 0;
isothermDenominator =0;

for jj = 1:omega
    isothermNumerator = isothermNumerator + (b.*P).^jj./(factorial(jj-1)).* ...
        ((1-jj*beta./vc)./(1-beta./vc)).^jj;
    isothermDenominator = isothermDenominator + (b.*P).^jj./(factorial(jj)).* ...
        ((1-jj*beta./vc)./(1-beta./vc)).^jj;
end
qa = isothermNumerator./(1+isothermDenominator);
end