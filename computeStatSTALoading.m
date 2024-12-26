function qa = computeStatSTALoading(P,T,b01,delU1,delU2,beta,kgate,cgate,sval,omega,vc)

b1 = b01.*exp(delU1./(8.314.*T));
b2 = b01.*exp(delU2./(8.314.*T));

isothermDenominator1 = 0;
isothermDenominator2 = 0;
isothermNumerator1   = 0;
isothermNumerator2   = 0;

for jj = 1:omega
% for jj = 1:round(vc./beta)
    isothermDenominator1 = isothermDenominator1 + ((b1.*P).^jj./(factorial(jj))).* ...
        ((1-jj*beta./vc)./(1-beta./vc)).^jj;
    isothermDenominator2 = isothermDenominator2 + ((b2.*P).^jj./(factorial(jj))).* ...
        ((1-jj*beta./vc)./(1-beta./vc)).^jj;

    isothermNumerator1 = isothermNumerator1 + (b1.*P).^(jj)./(factorial(jj-1)).* ...
        ((1-jj*beta./vc)./(1-beta./vc)).^jj;

    isothermNumerator2 = isothermNumerator2 + (b2.*P).^(jj)./(factorial(jj-1)).* ...
        ((1-jj*beta./vc)./(1-beta./vc)).^jj;
end

% q1 = isothermNumerator1./(1+isothermDenominator1);
% q2 = isothermNumerator2./(1+isothermDenominator2);

% ytrans = ((1+isothermDenominator2)./(1+isothermDenominator1)) .* exp(-(kgate - T.*cgate)./(8.314.*T));
% ytrans = (q2./q1) .* exp(-(kgate - T.*cgate)./(8.314.*T));
% ytrans = exp((kgate - T.*cgate)./(8.314.*T));
% ytrans = ((1+b2.*P)./(1+b1.*P)).^omega .* exp((kgate - T.*cgate)./(8.314.*T));
% ytrans = (1+b2.*P)./(1+b1.*P) .* exp(-(kgate - T.*cgate)./(8.314.*T));
% sigma = ytrans.^sval./(1+ytrans.^sval);
% sigma = (0.5 + 0.5.*tanh(-0.5.*(kgate - T.*cgate)./(8.314.*T)));




% qa = (1-sigma).*q1 + sigma.*q2;
% qa = (isothermNumerator1 + exp(-(kgate - T.*cgate)./(8.314.*T)).*isothermNumerator2)./ ...
%     ((1+isothermDenominator1) + exp(-(kgate - T.*cgate)./(8.314.*T)).*(1+isothermDenominator2));

qa = (isothermNumerator1 + exp(-(kgate./T - cgate)).*isothermNumerator2)./ ...
    ((1+isothermDenominator1) + exp(-(kgate./T - cgate)).*(1+isothermDenominator2));

% if isnan(any(qa))
%     qa(:) = 0;
% end


% plot(P,sigma)
% plot(P,qa)
end