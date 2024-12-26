function qa = computeStatSTALoading2(P,T,b01,b02,delU1,delU2,beta,kgate,cgate,delvc,omega,vc,varargin)

if length(varargin) ~= 0
    gamma = varargin{1};
else
    gamma = 1;
end

% delvc2 = delvc.*vc;

betavdW = vc./omega;

vc1 = vc - delvc.*vc;

omega1 = floor(vc1./betavdW);


b1 = b01.*exp(delU1./(8.314.*T));
b2 = b02.*exp(delU2./(8.314.*T));

isothermDenominator1 = 0;
isothermDenominator2 = 0;
isothermNumerator1   = 0;
isothermNumerator2   = 0;
% figure
jj = 2;

if omega1>0
    % for jj = 1:min(omega1,(vc1./beta))
    % for jj = 1:1000
    while jj <= min(floor(vc1./beta),omega1)
    % while jj <= vc1./beta 
        isothermNumerator1 = isothermNumerator1 + ((b1.*P).^gamma).^(jj)./(factorial(jj-1)).* ...
            ((1-jj*beta./vc1)./(1-beta./vc1)).^jj;
        isothermDenominator1 = isothermDenominator1 + (((b1.*P).^gamma).^jj./(factorial(jj))).* ...
            ((1-jj*beta./vc1)./(1-beta./vc1)).^jj;
        jj = jj+1;
        %
        % scatter(jj,((b1.*P(end)).^jj./(factorial(jj))).* ...
        %     ((1-jj*beta./vc1)./(1-beta./vc1)).^jj,'or')
        %
        % scatter(jj,isothermDenominator1(end),'ok')
    end

    % end
end

jj = 2;
% for jj = 1:min(omega,(vc./beta))
% for jj = 1:1000
while jj <= min(floor(vc./beta),omega)
% while jj <= vc./beta 
    isothermNumerator2 = isothermNumerator2 + ((b2.*P).^gamma).^(jj)./(factorial(jj-1)).* ...
        ((1-jj*beta./vc)./(1-beta./vc)).^jj;
    isothermDenominator2 = isothermDenominator2 + (((b2.*P).^gamma).^jj./(factorial(jj))).* ...
        ((1-jj*beta./vc)./(1-beta./vc)).^jj;
    jj = jj+1;
end
% end

% omega = min(omega,floor(vc./beta));
% omega1 = min(omega1,floor(vc1./beta));

%
% qaNum = (isothermNumerator1./(1+isothermDenominator1)).*(1+isothermDenominator1).^omega1 + ...
%         exp(-(kgate./T - cgate)).* ...
%         (isothermNumerator2./(1+isothermDenominator2)).*(1+isothermDenominator2).^omega;
%
% qaDen = (1 + isothermDenominator1).^omega1 + ...
%     exp(-(kgate./T - cgate)).* ...
%         (1 + isothermDenominator2).^omega;

qa1 = (((b1.*P).^gamma+isothermNumerator1)./(1+(b1.*P).^gamma+isothermDenominator1));
qa2 = (((b2.*P).^gamma+isothermNumerator2)./(1+(b2.*P).^gamma+isothermDenominator2));

yval = ((1+(b2.*P).^gamma+isothermDenominator2))./((1+(b1.*P).^gamma+isothermDenominator1)).*exp(-(kgate.*(8.314) - T.*cgate.*(8.314))./(T.*(8.314)));

sigmaval = yval./(1+yval);
qa = (1-sigmaval).*qa1 + sigmaval.*qa2;

% actFun = exp(-(kgate.*8.314 - T.*cgate.*8.314)./(8.314.*T));

% qa = (isothermNumerator1 + actFun.*isothermNumerator2)./((1+isothermDenominator1) + actFun.*(1+isothermDenominator2));
%


% subplot(1,2,1)
% hold on
% plot(P, (isothermNumerator1))
% plot(P, (isothermNumerator2).*exp(-(kgate./T - cgate)))
% set(gca,'YScale','log')
% subplot(1,2,2)
% hold on
% plot(P, (1+isothermDenominator1))
% plot(P, (1+isothermDenominator2).*exp(-(kgate./T - cgate)))
% set(gca,'YScale','log')

end