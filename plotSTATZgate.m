Pvals = linspace(0,2,1000);
qvals = zeros(1,length(Pvals));
Tvals = [273, 298, 308];
omega = 48;
beta = 6.0005e+01;
b01 = 4.5654e-04;
delU1 = 2.2271e+04;
vc1 = 2918;
vc2 = 3034;
gamma = 0.1;

kvals = 15;
cvals = 20;

figure(99)
for mm = 1:length(Tvals)
    T = Tvals(mm);
    for jj = 1:length(kvals)
        kgate = kvals(jj);
        for kk = 1:length(cvals)
            cgate = cvals(kk);
            for ii = 1:length(qvals)
                qvals(ii) = computeStatZGATELoading(Pvals(ii),T,b01,delU1,beta,omega,kgate,cgate,gamma,vc1,vc2);
            end
            subplot(2,2,1)
            plot(Pvals,qvals)
            hold on
            subplot(2,2,2)
            vcvals = (vc1+vc2)./2 + (vc2-vc1)./2.*tanh(kgate.*Pvals - (cgate*(exp(-10000./(8.314.*T)))));
            plot(Pvals,vcvals)
            hold on
            subplot(2,2,3)
            plot(Pvals,gradient(qvals', max(Pvals)./length(Pvals)))
            hold on
            subplot(2,2,4)
            scatter(T,-cgate*(exp(-10000./(8.314.*T))))
            hold on
        end
    end
end