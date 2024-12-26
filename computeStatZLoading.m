function qa = computeStatZLoading(P,T,b01,delU1,beta,omega,vc)

b = b01.*exp(delU1./(8.314.*T));

isothermNumerator = 0;
isothermDenominator = 0;

jj = 2;
% for jj = 1:omega
while jj <= min(vc./beta,omega)
    isothermNumerator = isothermNumerator + (b.*P).^jj./(factorial(jj-1)).* ...
        ((1-jj*beta./vc)./(1-beta./vc)).^jj;
    isothermDenominator = isothermDenominator + (b.*P).^jj./(factorial(jj)).* ...
        ((1-jj*beta./vc)./(1-beta./vc)).^jj;
    jj = jj+1;
end

qa = (b.*P+isothermNumerator)./(1+b.*P+isothermDenominator);
end