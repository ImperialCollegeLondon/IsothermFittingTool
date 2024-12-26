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
% parameters for gate opening model
%
% Last modified:
% - 2023-06-11, HA: Initial creation
%
% Input arguments:
%
% Output arguments:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
% close all;
fprintf('Load isotherm file for species A from fittingResults');
uiopen;
isothermDataA = isothermData;
omegaA = isothermDataA.isothermParameters(1,1);
betaA =  isothermDataA.isothermParameters(2,1);
b01A =   isothermDataA.isothermParameters(3,1);
delU1A = isothermDataA.isothermParameters(4,1);
delU2A = isothermDataA.isothermParameters(5,1);
kgateA = isothermDataA.isothermParameters(6,1);
cgateA = isothermDataA.isothermParameters(7,1);

fprintf('Load isotherm file for species B from fittingResults');
uiopen;
isothermDataB = isothermData;
omegaB = isothermDataB.isothermParameters(1,1);
betaB =  isothermDataB.isothermParameters(2,1);
b01B =   isothermDataB.isothermParameters(3,1);
delU1B = isothermDataB.isothermParameters(4,1);
if length(isothermDataB.isothermParameters(:,1)) > 4
    delU2B = isothermDataB.isothermParameters(5,1);
    kgateB = isothermDataB.isothermParameters(6,1);
    cgateB = isothermDataB.isothermParameters(7,1);
else
    delU2B = delU1B;
    kgateB = 1;
    cgateB = 1;
end
LitData = [];
vc = isothermData.CageVolume;
vm =  isothermData.MicroporeVolume;
Na = 6.022e20; % Avogadros constant [molecules/mmol];

%%
P = linspace(0,max(isothermDataA.experiment(:,1))*1.1,10000);
T = 303.15;
% yA = 0.66;
% yA = 0.76;
% yA = 0.83;
yA = 0.5;
[qA, qB, qT]  = computeStatZGOLoadingBinary(P,T,b01A,delU1A,delU2A,kgateA,cgateA,betaA,omegaA,b01B,delU1B,delU2B,kgateB,cgateB,betaB,omegaB,vc,yA);
qAmolkg = qA.*vm./(vc.*Na);
qBmolkg = qB.*vm./(vc.*Na);
qTmolkg = qAmolkg+qBmolkg;
% [qA, qB, qT]  = computeStatZLoadingBinary(P,T,b01A,delU1A,betaA,omegaA,b01A,delU1A,betaA,omegaA,vc,yA);
% [qA, qB, qT]  = computeStatZLoadingBinary(P,T,6*0.75e-4,5.1*4200,64.5,12,3.5*0.75e-4,5*4200,77,10,776,yA);
% [qA, qB, qT] = computeStatZLoadingBinary(P,T,6e-7,5.1*4184,64.5,12,3.5e-7,5*4184,77,10,776,yA);
qc = computeStatZGATELoading2(P,T,b01A,delU1A,delU2A,betaA,kgateA,cgateA,omegaA,vc);
qd = computeStatZGATELoading2(P,T,b01B,delU1B,delU2B,betaB,kgateB,cgateB,omegaB,vc);
% % qc = computeStatZLoading(linspace(0,P,length(yA)),T,6e-7,5.1*4184,64.5,12,776); qd = computeStatZLoading(linspace(0,P,length(yA)),T,3.5e-7,5*4184,77,10,776);
% qd = qc;


% LitData = [206.0,350.4,505.0,656.7,804.5,904,1004;
%     0.27,0.35,0.48,0.60,0.94,1.40,2.02;
%     0.02,0.04,0.05,0.09,0.06,0.10,0.19]; % 0.66
% 
% LitData = [203.8,348.0,503.0,654.6,802.5,901.8,1001.5;
%     0.28,0.40,0.52,0.70,1.36,2.03,2.32;
%     0.02,0.03,0.04,0.06,0.09,0.07,0.15]; % 0.76
% % %
% LitData = [205.3, 349.6, 504.4,655.7, 803.4, 903.1;
%  0.30, 0.43,0.58,0.86,1.93,2.24;
%  0.01,0.01,0.00,0.00,0.00,0.00]; % 0.83


figure
hold on
plot(P,qAmolkg,'r','MarkerFaceColor','none','MarkerEdgeColor','none','LineWidth',2,'DisplayName','Mixed Gas A')
plot(P,qBmolkg,'b','MarkerFaceColor','none','MarkerEdgeColor','none','LineWidth',2,'DisplayName','Mixed Gas B')
plot(P,qTmolkg,'--k','MarkerFaceColor','none','MarkerEdgeColor','none','LineWidth',2,'DisplayName','total')
% legend('Location','northwest')
if ~isempty(LitData)
    scatter(LitData(1,:)./100,LitData(2,:),70,'r','filled','DisplayName','Exp A')
    scatter(LitData(1,:)./100,LitData(3,:),70,'b','filled','DisplayName','Exp B')
end
xlabel('Pressure [bar]');
ylabel('Amount adsorbed [mol/kg]');
xlim([0.1 P(end)])
box on
set(gca,'YScale','linear','XScale','log','FontSize',15,'LineWidth',1)
grid on;
% axis square
set(gcf,'units','inch','position',[7,5,5,5])
%
figure
hold on
plot(P,qc.*vm./(vc.*Na),'r','MarkerFaceColor','none','MarkerEdgeColor','none','LineWidth',2,'DisplayName','Pure Gas A')
plot(isothermDataA.experiment(find(isothermDataA.experiment(:,3)==T),1),isothermDataA.experiment(find(isothermDataA.experiment(:,3)==T),2).*vm./(vc.*Na),'or','MarkerFaceColor','r','MarkerEdgeColor','r','HandleVisibility','off','MarkerSize',8)
plot(P,qd.*vm./(vc.*Na),'b','MarkerFaceColor','none','MarkerEdgeColor','none','LineWidth',2,'DisplayName','Pure Gas B')
plot(isothermDataB.experiment(find(isothermDataB.experiment(:,3)==T),1),isothermDataB.experiment(find(isothermDataB.experiment(:,3)==T),2).*vm./(vc.*Na),'ob','MarkerFaceColor','b','MarkerEdgeColor','b','HandleVisibility','off','MarkerSize',8)
legend('Location','northwest')
xlabel('Pressure [bar]');
ylabel('Amount adsorbed [mol/kg]');
xlim([0.1 P(end)])
box on
set(gca,'YScale','linear','XScale','log','FontSize',15,'LineWidth',1)
grid on;
% axis square
set(gcf,'units','inch','position',[0,5,5,5])
%
figure
delUA = (delU1A+delU2A)./2 + (delU2A-delU1A)./2.*tanh((kgateA.*((b01A.*exp(delU1A./(8.314.*T)))./(1+b01A.*exp(delU1A./(8.314.*T)).*P)).*P - cgateA)./2);
plot(P,delUA,'r','MarkerFaceColor','none','MarkerEdgeColor','none','LineWidth',2,'DisplayName','Pure Gas A')
hold on
delUB = (delU1B+delU2B)./2 + (delU2B-delU1B)./2.*tanh((kgateB.*((b01B.*exp(delU1B./(8.314.*T)))./(1+b01B.*exp(delU1B./(8.314.*T)).*P)).*P - cgateB)./2);
plot(P,delUB,'b','MarkerFaceColor','none','MarkerEdgeColor','none','LineWidth',2,'DisplayName','Pure Gas B')
xlabel('P [bar]')
ylabel('deltaU [J/mol]')
hold on
box on
set(gca,'YScale','linear','XScale','log','FontSize',15,'LineWidth',1)
grid on;
xlim([0.1 P(end)])
set(gcf,'units','inch','position',[10,5,5,5])
% legend('Location','northwest')

