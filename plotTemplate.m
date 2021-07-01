figure
load('C:\Users\azxan\Documents\GitHub\IsothermFittingTool\isothermData\AC_Quanta_1_ads.mat')
scatter(fitData(:,1),fitData(:,2),10,'MarkerEdgeColor','#E5383B')
hold on
load('C:\Users\azxan\Documents\GitHub\IsothermFittingTool\isothermData\AC_Quanta_2_ads.mat')
scatter(fitData(:,1),fitData(:,2),10,'MarkerEdgeColor','#008FF5')
hold on
load('C:\Users\azxan\Documents\GitHub\IsothermFittingTool\isothermData\AC_Ronny_all.mat')
scatter(fitData(:,1),fitData(:,2),10,'MarkerEdgeColor','#6C757D')
hold on
load('C:\Users\azxan\Documents\GitHub\IsothermFittingTool\isothermData\AC_Ronny.mat')
scatter(fitData(:,1),fitData(:,2),10,'MarkerEdgeColor','g')
box on; grid on; axis square
set(gca,'YScale','linear','XScale','linear','FontSize',10,'LineWidth',1)
set(gcf,'units','inch','position',[0,0,3.3,3])
xlabel('P [bar]');
ylabel('q* [mol/kg]')

legend('AC sample 1','AC sample 2','AC literature','box','off','location','northwest')
