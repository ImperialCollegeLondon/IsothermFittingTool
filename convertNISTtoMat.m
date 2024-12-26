% selpath = uigetdir;
% addpath(genpath(selpath));
% cd(selpath);
% dirData = dir('**/*.csv');
% [path FolderName ext] = fileparts(dirData(1).folder);
fitDataAll = [];
finalLoading = 1e3;
% for ii = 1:size(dirData,1)
%     tempData = readmatrix([dirData(ii).name]);
%     Tval = tempData(find(~isnan(tempData(:,2)), 1, 'first'),2);
%     initInd = find(~isnan(tempData(:,3)), 1, 'first');
%     Pvals = tempData(initInd:end,2);
%     qvals = tempData(initInd:end,4);
%     Tvals = ones(length(qvals),1).*Tval;
% 
%     fitDataAll = [fitDataAll; Pvals qvals Tvals];
% end

fitDataAll = sortrows(fitDataAll,2);
fitDataAll = sortrows(fitDataAll,3);

x = fitDataAll(:,1);
z = fitDataAll(:,2);
y = fitDataAll(:,3);

temperatureValues = unique(y);
qRefIndexTemp = zeros(length(temperatureValues),1);
for ii = 1:length(temperatureValues)
    qRefIndexTemp(ii,1) = find(y == temperatureValues(ii),1,'first');
    qRefIndexTemp(ii,2) = find(y == temperatureValues(ii),1,'last');
    qRefIndexTemp(ii,3) = qRefIndexTemp(ii,1) - 1 + find(z(qRefIndexTemp(ii,1):qRefIndexTemp(ii,2)) < finalLoading,1,'last');
    % x = x
end

desFlag = [];
fitDataDes = [];
fitData = [];

for mm = 1:length(temperatureValues)
    maxVal = x(qRefIndexTemp(mm,1));
    for jj = qRefIndexTemp(mm,1)+1:qRefIndexTemp(mm,3)
        if x(jj,1) >= maxVal && x(jj,1) >= x(jj-1,1)
            maxVal = x(jj,1);
            desFlag(jj) = false;
            fitData = [fitData; fitDataAll(jj,:)];
        else
            desFlag(jj) = true;
            fitDataDes = [fitDataDes; fitDataAll(jj,:)];
        end
    end
end

desFlag = logical(desFlag);


% figure
subplot(2,3,1)
scatter(fitData(:,1),fitData(:,2),'ob')
hold on
subplot(2,3,2)
if ~isempty(fitDataDes)
scatter(fitDataDes(:,1),fitDataDes(:,2),'or')
end
hold on
subplot(2,3,3)
hold on
if ~isempty(fitDataDes)
scatter(fitDataDes(:,1),fitDataDes(:,2),'or')
end
scatter(fitData(:,1),fitData(:,2),'ob')
subplot(2,3,4)
scatter(fitData(:,1),fitData(:,2),'ob')
hold on
set(gca,'XScale','log')
subplot(2,3,5)
if ~isempty(fitDataDes)
scatter(fitDataDes(:,1),fitDataDes(:,2),'or')
end
hold on
set(gca,'XScale','log')
subplot(2,3,6)
hold on
if ~isempty(fitDataDes)
scatter(fitDataDes(:,1),fitDataDes(:,2),'or')
end
scatter(fitData(:,1),fitData(:,2),'ob')
set(gca,'XScale','log')

figure
scatter(fitDataAll(qRefIndexTemp(4,1):qRefIndexTemp(4,3),1),fitDataAll(qRefIndexTemp(4,1):qRefIndexTemp(4,3),2),'ob')
set(gca,'XScale','log','YScale','log')
% fitDataDes = sortrows(fitDataDes,1);
% fitDataDes = sortrows(fitDataDes,3);
% Outs = isoutlier(fitDataDes(:,2));
% scatter(fitDataDes(Outs,1),fitDataDes(Outs,2),'yellow','filled')


cd('/Users/ha3215/Documents/GitHub/ERASE/IsothermFittingTool');
save(['..',filesep,'IsothermFittingTool',filesep','isothermData',filesep,FolderName,'_auto_',char(datetime('now','Format','ddMMyy_HHmm'))],"fitData","fitDataDes")


