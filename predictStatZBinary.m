%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Imperial College London, United Kingdom
% Multifunctional Nanomaterials Laboratory / Complex Porous Media
% Laboratory
%
% Project:  PhD
% Year:     2023
% MATLAB:   R2020a
% Authors:  Hassan Azzan (HA)
%
% Purpose:
% Predicts binary adsorption isotherms using statistical isotherm model
% parameters for zeolites
%
% Last modified:
% - 2023-06-11, HA: Initial creation
%
% Input arguments:
%
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all;
fprintf('Load isotherm file for species A from fittingResults');
uiopen;
isothermDataA = isothermData;
omegaA = isothermDataA.isothermParameters(1,1);
betaA = isothermDataA.isothermParameters(2,1);
b01A = isothermDataA.isothermParameters(3,1);
delU1A = isothermDataA.isothermParameters(4,1);

fprintf('Load isotherm file for species B from fittingResults');
uiopen;
isothermDataB = isothermData;
omegaB = isothermDataB.isothermParameters(1,1);
betaB = isothermDataB.isothermParameters(2,1);
b01B = isothermDataB.isothermParameters(3,1);
delU1B = isothermDataB.isothermParameters(4,1);
vc = isothermData.CageVolume;
vm =  isothermData.MicroporeVolume;
Na = 6.022e20; % Avogadros constant [molecules/mmol];

P = 3.2;
T = 298.15;
yA = linspace(0,1,500);
% yA = 0.25;
% yA = 0.15;
[qA, qB, qT]  = computeStatZLoadingBinary(P,T,b01A,delU1A,betaA,omegaA,b01B,delU1B,betaB,omegaB,vc,yA);
qAmolkg = qA.*vm./(vc.*Na);
qBmolkg = qB.*vm./(vc.*Na);
qTmolkg = qAmolkg+qBmolkg;
% [qA, qB, qT]  = computeStatZLoadingBinary(P,T,b01A,delU1A,betaA,omegaA,b01A,delU1A,betaA,omegaA,vc,yA);
% [qA, qB, qT]  = computeStatZLoadingBinary(P,T,6*0.75e-4,5.1*4200,64.5,12,3.5*0.75e-4,5*4200,77,10,776,yA);
% [qA, qB, qT] = computeStatZLoadingBinary(P,T,6e-7,5.1*4184,64.5,12,3.5e-7,5*4184,77,10,776,yA);
qc = computeStatZLoading(linspace(0,P,length(yA)),T,b01A,delU1A,betaA,omegaA,vc); qd = computeStatZLoading(linspace(0,P,length(yA)),T,b01B,delU1B,betaB,omegaB,vc);
% % qc = computeStatZLoading(linspace(0,P,length(yA)),T,6e-7,5.1*4184,64.5,12,776); qd = computeStatZLoading(linspace(0,P,length(yA)),T,3.5e-7,5*4184,77,10,776);
% qd = qc;

figure
plot(yA,qAmolkg,'r','MarkerFaceColor','none','MarkerEdgeColor','none','LineWidth',2,'DisplayName','Gas A')
hold on
plot(yA,qBmolkg,'b','MarkerFaceColor','none','MarkerEdgeColor','none','LineWidth',2,'DisplayName','Gas A')
plot(yA,qTmolkg,'k','MarkerFaceColor','none','MarkerEdgeColor','none','LineWidth',2,'DisplayName','Gas A')
xlabel('$$y_{A}$$ [-]');
ylabel('Amount adsorbed [mol/kg]');
xlim([0 1])
box on
legend('Location','best')
set(gca,'YScale','linear','XScale','linear','FontSize',15,'LineWidth',1)
grid on; axis square
set(gcf,'units','inch','position',[0,5,5,5])
% 
% figure
% scatter(linspace(0,P,length(yA)),qc,5,'filled','MarkerFaceColor','r','MarkerFaceAlpha',0.1,'MarkerEdgeColor','r','LineWidth',0.8,'DisplayName','Gas A')
% hold on
% scatter(linspace(0,P,length(yA)),qd,6,'filled','MarkerFaceColor','b','MarkerFaceAlpha',0.1,'MarkerEdgeColor','b','LineWidth',0.8,'DisplayName','Gas B')
% xlabel('P [bar]');
% ylabel('Amount adsorbed [molecules/supercage]');
% xlim([0 P])
% box on
% legend('Location','best')
% set(gca,'YScale','linear','XScale','linear','FontSize',15,'LineWidth',1)
% grid on; axis square
% set(gcf,'units','inch','position',[5,5,5,5])
% 
% figure
% scatter(yA,qA./(qT),5,'filled','MarkerFaceColor','r','MarkerFaceAlpha',0.1,'MarkerEdgeColor','r','LineWidth',0.8)
% hold on
% plot([0 1],[0 1],'k','LineStyle','--','LineWidth',1)
% xlabel('$$y_{A}$$ [-]');
% ylabel('$$x_{A}$$ [-]');
% xlim([0 1])
% box on
% % legend('Location','best')
% set(gca,'YScale','linear','XScale','linear','FontSize',15,'LineWidth',1)
% grid on; axis square
% set(gcf,'units','inch','position',[10,0,5,5])
% 
% figure
% selectivityA = qA./qB.*(1-yA')./yA';
% scatter(yA,selectivityA,5,'filled','MarkerFaceColor','r','MarkerFaceAlpha',0.1,'MarkerEdgeColor','r','LineWidth',0.8)
% hold on
% ylabel('$$\alpha_{A}$$ [-]');
% xlabel('$$y_{A}$$ [-]');
% xlim([0 1])
% box on
% % legend('Location','best')
% set(gca,'YScale','linear','XScale','linear','FontSize',15,'LineWidth',1)
% grid on; axis square
% set(gcf,'units','inch','position',[10,0,5,5])
% 
% figure
% [qADES, qBDES, qTDES]  = computeStatZLoadingBinary(1,T,b01A,delU1A,betaA,omegaA,b01B,delU1B,betaB,omegaB,vc,yA);
% WCA = qA-qADES;
% WCB = qB-qBDES;
% yangFOM = selectivityA.*(WCA./WCB);
% scatter(yA,yangFOM,5,'filled','MarkerFaceColor','r','MarkerFaceAlpha',0.1,'MarkerEdgeColor','r','LineWidth',0.8)
% hold on
% ylabel('YANG [-]');
% xlabel('$$y_{A}$$ [-]');
% xlim([0 1])
% box on
% % legend('Location','best')
% set(gca,'YScale','linear','XScale','linear','FontSize',15,'LineWidth',1)
% grid on; axis square
% set(gcf,'units','inch','position',[10,0,5,5])