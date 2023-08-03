function q  = computeUNIV6EDF(epsilon, a1, a2, a3, e01, e02, e03, e04, m1, m2, m3, m4, hfg)

a4 = 1-a1-a2-a3;

E1 = exp(e01)-1;
E2 = exp(e02)-1;
E3 = exp(e03)-1;
E4 = exp(e04)-1;
% hfg = 43.988*1000; % latent heat of vaporization of water at 298 K [J/mol]
% hfg = 8.5*1000; % latent heat of vaporization of methane [J/mol]
% hfg = 0;



term1 = a1.*(exp((epsilon-hfg-E1)./(m1))./(m1.*((1+exp((epsilon-hfg-E1)./(m1))).^2)));
term2 = a2.*(exp((epsilon-hfg-E2)./(m2))./(m2.*((1+exp((epsilon-hfg-E2)./(m2))).^2)));
term3 = a3.*(exp((epsilon-hfg-E3)./(m3))./(m3.*((1+exp((epsilon-hfg-E3)./(m3))).^2)));
if e04 == 0
    term4 = zeros(1,length(epsilon));
else
    term4 = a4.*(exp((epsilon-hfg-E4)./(m4))./(m4.*((1+exp((epsilon-hfg-E4)./(m4))).^2)));
end

q = [term1;term2;term3;term4];

end