function qa = computeStatZELoading(P,T,b01,delU1,beta,omega,ebyk,vc)

R = 8.314;

% b = b01.*exp((delU1+ebyk.*beta./vc)./(8.314.*T));
b = b01.*exp((delU1)./(8.314.*T));

isothermNumerator = 0;
isothermDenominator = 0;
for jj = 1:omega
    if jj*beta < vc
        % isothermNumerator = isothermNumerator + (b.*P).^jj./(factorial(jj-1)).* ...
        %     (((1-jj*beta./vc)./(1-beta./vc)).^jj).*exp(jj.*beta.*ebyk./(vc.*T));
        % isothermDenominator = isothermDenominator + (b.*P).^jj./(factorial(jj)).* ...
        %     (((1-jj*beta./vc)./(1-beta./vc)).^jj).*exp(jj.*beta.*ebyk./(vc.*T));

        isothermNumerator = isothermNumerator + ((b.*P).^(R.*T./ebyk)).^(jj)./(factorial(jj-1)).* ...
            ((1-jj*beta./vc)./(1-beta./vc)).^jj;
        isothermDenominator = isothermDenominator + ((b.*P).^(R.*T./ebyk)).^(jj)./(factorial(jj)).* ...
            ((1-jj*beta./vc)./(1-beta./vc)).^jj;
    else
        break
    end
end
qa = (isothermNumerator)./(1+isothermDenominator);

end