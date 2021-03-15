function [outScatter]=generateUncertaintySpread(x,y,isothermModel,parVals,conRange95)
clc
nPoints = 20;
Pvals=linspace(0,max(x),1000);
Tvals=unique(y);

switch isothermModel
    case 'DSL'
        qs1 = parVals(1);
        qs2 = parVals(2);
        b01 =  parVals(3);
        b02 =  parVals(4);
        delU1 =  parVals(5);
        delU2= parVals(6);
        
        qeqUnc=zeros(3,nPoints*length(Pvals),length(Tvals));
        for mm = 1:length(Tvals)
            for jj = 1:length(Pvals)
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj+1))
                    qeqUnc(1,kk,mm)=Pvals(jj);
                    qeqUnc(2,kk,mm)=Tvals(mm);
                end
            end
        end
        
        lhsMat = 2*lhsdesign(nPoints,6)-1;
        
        for mm = 1:length(Tvals)
            for jj = 1:length(Pvals)
                hh = 0;
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj))
                    hh=hh+1;
                    qs1_unc = qs1+conRange95(1)*lhsMat(hh,1);
                    qs2_unc = qs2+conRange95(2)*lhsMat(hh,2);
                    b01_unc = b01+conRange95(3)*lhsMat(hh,3);
                    b02_unc = b02+conRange95(4)*lhsMat(hh,4);
                    delU1_unc = delU1+conRange95(5)*lhsMat(hh,5);
                    delU2_unc = delU2+conRange95(6)*lhsMat(hh,6);
                    P = qeqUnc(1,kk,mm);
                    T = qeqUnc(2,kk,mm);
                    qeqUnc(3,kk,mm)=qs1_unc.*(b01_unc.*P.*exp(delU1_unc./(8.314.*T)))./(1+(b01_unc.*P.*exp(delU1_unc./(8.314.*T)))) ...
                        + qs2_unc.*(b02_unc.*P.*exp(delU2_unc./(8.314.*T)))./(1+(b02_unc.*P.*exp(delU2_unc./(8.314.*T))));
                end
            end
        end
        
        outScatter=[qeqUnc(1,:);qeqUnc(2,:); qeqUnc(3,:)];
        outScatter=outScatter';
        
    case 'DSS'
        qs1 = parVals(1);
        qs2 = parVals(2);
        b01 =  parVals(3);
        b02 =  parVals(4);
        delU1 =  parVals(5);
        delU2= parVals(6);
        gamma= parVals(7);
        
        qeqUnc=zeros(3,nPoints*length(Pvals),length(Tvals));
        for mm = 1:length(Tvals)
            for jj = 1:length(Pvals)
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj+1))
                    qeqUnc(1,kk,mm)=Pvals(jj);
                    qeqUnc(2,kk,mm)=Tvals(mm);
                end
            end
        end
        
        lhsMat = 2*lhsdesign(nPoints,7)-1;
        
        for mm = 1:length(Tvals)
            for jj = 1:length(Pvals)
                hh = 0;
                for kk = (nPoints*(jj-1)+1):(nPoints*(jj))
                    hh=hh+1;
                    qs1_unc = qs1+conRange95(1)*lhsMat(hh,1);
                    qs2_unc = qs2+conRange95(2)*lhsMat(hh,2);
                    b01_unc = b01+conRange95(3)*lhsMat(hh,3);
                    b02_unc = b02+conRange95(4)*lhsMat(hh,4);
                    delU1_unc = delU1+conRange95(5)*lhsMat(hh,5);
                    delU2_unc = delU2+conRange95(6)*lhsMat(hh,6);
                    gamma_unc = gamma+conRange95(7)*lhsMat(hh,7);
                    P = qeqUnc(1,kk,mm);
                    T = qeqUnc(2,kk,mm);
                    qeqUnc(3,kk,mm)=qs1_unc.*(b01_unc.*P.*exp(delU1_unc./(8.314.*T))).^gamma_unc./(1+(b01_unc.*P.*exp(delU1_unc./(8.314.*T))).^gamma_unc) ...
                        + qs2_unc.*(b02_unc.*P.*exp(delU2_unc./(8.314.*T))).^gamma_unc./(1+(b02_unc.*P.*exp(delU2_unc./(8.314.*T))).^gamma_unc);
                end
            end
        end
        
        outScatter=[qeqUnc(1,:);qeqUnc(2,:); qeqUnc(3,:)];
        outScatter=outScatter';
end
end