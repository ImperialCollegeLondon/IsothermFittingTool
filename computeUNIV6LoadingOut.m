function q = computeUNIV6LoadingOut(p,T, qs, a1, a2, a3, e01, e02, e03, e04, m1, m2, m3, m4, ps)

a4 = 1-a1-a2-a3;
R = 8.314;
x = p./ps;

term1 = a1.*(x.*exp(e01./(R.*T))).^(R.*T./m1)./(1+(x.*exp(e01./(R.*T))).^(R.*T./m1));
term2 = a2.*(x.*exp(e02./(R.*T))).^(R.*T./m2)./(1+(x.*exp(e02./(R.*T))).^(R.*T./m2));
term3 = a3.*(x.*exp(e03./(R.*T))).^(R.*T./m3)./(1+(x.*exp(e03./(R.*T))).^(R.*T./m3));
term4 = a4.*(x.*exp(e04./(R.*T))).^(R.*T./m4)./(1+(x.*exp(e04./(R.*T))).^(R.*T./m4));

if a4 > 0
    q = qs.*(term1+term2+term3+term4);
else
    q = 0;
end

end