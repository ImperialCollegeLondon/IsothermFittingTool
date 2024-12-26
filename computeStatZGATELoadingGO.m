function qa = computeStatZGATELoadingGO(P,T,b01,b02,delU1,delU2,beta,kgate,cgate,delvc,omega,vc)



% delvc = delvc.*vc;

betavdW = vc./omega;

vc1 = vc - delvc.*vc;

omega1 = ceil(vc1./betavdW);

% sigma = (1)./2 + (1)./2.*tanh((kgate.*((b01.*P.*exp(delU1./(8.314.*T)))./(1+b01.*exp(delU1./(8.314.*T)).*P)) - cgate)./2);
sigma = (1)./2 + (1)./2.*tanh((kgate.*((b01.*P.*exp(delU1./(8.314.*T))).*P) - cgate)./2);

b1 = b01.*exp(delU1./(8.314.*T));
b2 = b02.*exp(delU2./(8.314.*T));


isothermNumerator1 = 0;
isothermDenominator1 = 0;
isothermNumerator2 = 0;
isothermDenominator2 = 0;


if omega1>0
    % for jj = 1:min(omega1,vc1./beta)
    for jj = 1:omega1
        if (jj*beta <= vc1)
            isothermNumerator1 = isothermNumerator1 + (b1.*P).^(jj)./(factorial(jj-1)).* ...
                ((1-jj*beta./vc1)./(1-beta./vc1)).^jj;
            isothermDenominator1 = isothermDenominator1 + ((b1.*P).^jj./(factorial(jj))).* ...
                ((1-jj*beta./vc1)./(1-beta./vc1)).^jj;
        end

    end
end

% for jj = 1:min(omega,vc./beta)
for jj = 1:omega
    if (jj*beta <= vc)
        isothermNumerator2 = isothermNumerator2 + (b2.*P).^(jj)./(factorial(jj-1)).* ...
            ((1-jj*beta./vc)./(1-beta./vc)).^jj;
        isothermDenominator2 = isothermDenominator2 + ((b2.*P).^jj./(factorial(jj))).* ...
            ((1-jj*beta./vc)./(1-beta./vc)).^jj;
    end
end


qa1 = isothermNumerator1./(1+isothermDenominator1);
qa2 = isothermNumerator2./(1+isothermDenominator2);
qa = (1-sigma).*qa1 + sigma.*qa2;

if isnan(qa(:))
    qa = zeros(length(P),1);
end

end