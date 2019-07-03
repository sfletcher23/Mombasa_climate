%% Optimal policies plots: Fig 4

load('/Users/sarahfletcher/Dropbox (MIT)/Mombasa_Climate/SDP/results67847_19_Mar_2018_21_47_26_base.mat')
totalCostTime = totalCostTime{1};
% Regret for last time period
bestOption = 0;

P_regret = [68 78 88];
totalCost = squeeze(sum(totalCostTime(:,:,1:3), 2));
meanCostPnow = zeros(length(P_regret),3);
for i = 1:length(P_regret)
    % Find simulaitons with this level of precip
    ind_P = P_state(:,end) == P_regret(i);
    % Get average cost of each infra option in that P level
    totalCostPnow = totalCost(ind_P,:);
    meanCostPnow(i,:) = mean(totalCostPnow,1);
end

bestInfraCost = min(meanCostPnow,[],2);
regret = meanCostPnow - repmat(bestInfraCost, 1,3);

f = figure;
font_size = 12;
bar([meanCostPnow; regret]/1E6)
hold on
line([3.5 3.5], [0 180],'Color', 'k')
xlim([.5 6.5])
xticklabels({'68', '78', '88', '68', '78', '88'})
yl = ylabel('M$')
yl.Position = yl.Position - [ .3 0 0];
xl = xlabel('P in 2090 [mm/month]')
xl.Position = xl.Position - [ 0 4 0];
title('Cost and Regret for Infrastructure Alternatives by 2090 P')
l = legend('Flexible', 'Large', 'Small')
l.Position = l.Position + [-.1 -.1 0 0.1]
legend('boxoff')
% FONT
allaxes = findall(f, 'type', 'axes');
set(allaxes,'FontSize', font_size)
set(findall(allaxes,'type','text'),'FontSize', font_size)
printsetup(f, 7, 5, 12, 1, 300, 1, 1, 'regret' )

