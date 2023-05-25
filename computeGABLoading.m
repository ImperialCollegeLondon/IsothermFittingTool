function qa = computeGABLoading(P,T,qs1,parC,parD,parF,parG)
E1 = parC - exp(parD.*T);
E29 = parF + parG.*T;
E10 = -44.38.*T + 57220;

c = exp((E1-E10)./(8.314.*T));
k = exp((E29-E10)./(8.314.*T));

qa = qs1.*k.*c.*P./((1-k.*P).*(1+(c-1).*k.*P));
end