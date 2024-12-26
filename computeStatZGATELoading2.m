function qa = computeStatZGATELoading2(P,T,b01,delU1,delU2,beta,kgate,cgate,omega,vc)

% delU = (delU1+delU2)./2 + (delU2-delU1)./2.*tanh(kgate.*b01.*exp(delU1./(8.314.*T)).*P - cgate);
delU = (delU1+delU2)./2 + (delU2-delU1)./2.*tanh((kgate.*((b01.*exp(delU1./(8.314.*T)))./(1+b01.*exp(delU1./(8.314.*T)).*P)).*P - cgate)./2);
% delU = (delU1+delU2)./2 + (delU2-delU1)./2.*tanh((kgate.*((b01.*exp(delU1./(8.314.*T)))).*P - cgate)./2);
% delU = (delU1+delU2)./2 + (delU2-delU1)./2.*tanh((kgate.*(b01.*exp(delU1./(8.314.*T)))./(1+b01.*exp(delU1./(8.314.*T)).*P).*P - cgate./(kgate.*(b01.*exp(delU1./(8.314.*T)))./(1+b01.*exp(delU1./(8.314.*T)).*P)));
% k = kgate.*(b01.*exp(delU1./(8.314.*T)).*P)./(1+b01.*exp(delU1./(8.314.*T)).*P);
% P0 = cgate;
% delU = ((delU1+delU2)./2 + (delU2-delU1)./2)./(1+exp(-k.*(P-P0)));

b = b01.*exp(delU./(8.314.*T));
% vc = (vc1+vc2)./2 + (vc2-vc1)./2.*tanh(kgate.*P - cgate*(exp(-10000./(8.314.*T))));

isothermNumerator = 0;
isothermDenominator = 0;

% omega = vc./beta;

for jj = 1:omega
    isothermNumerator = isothermNumerator + ((b.*P).^jj./(factorial(jj-1))).* ...
        ((1-jj*beta./vc)./(1-beta./vc)).^jj;
    isothermDenominator = isothermDenominator + ((b.*P).^jj./(factorial(jj))).* ...
        ((1-jj*beta./vc)./(1-beta./vc)).^jj;
end
qa = isothermNumerator./(1+isothermDenominator);
% qa = (isothermNumerator./(1+isothermDenominator)).^(1./gamma);

end