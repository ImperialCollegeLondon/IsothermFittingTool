function qa = computeStatZLoading(P,T,b01,delU1,beta,omega,vc)

b = b01.*exp(delU1./(8.314.*T));

isothermNumerator = 0;
isothermDenominator = 0;

for jj = 1:omega+1
    if jj*beta < vc
        isothermNumerator = isothermNumerator + (b.*P).^jj./(factorial(jj-1)).* ...
            ((1-jj*beta./vc)).^jj;
        isothermDenominator = isothermDenominator + (b.*P).^jj./(factorial(jj)).* ...
            ((1-jj*beta./vc)).^jj;
    else
        break
    end
end
qa = isothermNumerator./(1+isothermDenominator);
end