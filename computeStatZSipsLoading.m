function qa = computeStatZSipsLoading(P,T,b01,delU1,beta,omega, gamma, vc)

b = b01.*exp(delU1./(8.314.*T));

isothermNumerator = 0;
isothermDenominator = 0;
for jj = 2:omega
    if jj*beta < vc
        isothermNumerator = isothermNumerator + ((b.*P).^gamma).^jj./(factorial(jj-1)).* ...
            ((1-jj*beta./vc)).^jj;
        isothermDenominator = isothermDenominator + ((b.*P).^gamma).^jj./(factorial(jj)).* ...
            ((1-jj*beta./vc)).^jj;
    else
        break
    end
end

qa = ((b.*P).^gamma+isothermNumerator)./(1+(b.*P).^gamma+isothermDenominator);


% kb = 1.38e-23;
% 
% for jj = 2:omega
%     if jj*beta < vc
%         isothermNumerator = isothermNumerator + (b.*P).^jj./(factorial(jj-1)).* ...
%             ((1-jj*beta./vc)).^jj.*exp((jj.*beta.*gamma)./(vc.*8.314.*T));
%         isothermDenominator = isothermDenominator + (b.*P).^jj./(factorial(jj)).* ...
%             ((1-jj*beta./vc)).^jj.*exp((jj.*beta.*gamma)./(vc.*8.314.*T));
%     else
%         break
%     end
% end
% 
% qa = (b.*P+isothermNumerator)./(1+b.*P+isothermDenominator);


end