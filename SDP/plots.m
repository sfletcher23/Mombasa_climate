%% Plots plots plots !!!

decade = {'2001-2020', '2021-2040', '2041-2060', '2061-2080', '2081-2100'};
%% Analyze runoff and shortage
if false
    
avgMAR = cellfun(@(x) mean(mean(x,2)), runoff);
avgStd = cellfun(@(x) mean(std(x,0,2)), runoff);


% Boxplots by P
avgMarNow = avgMAR(:,1:28,5);
avgStdNow = avgStd(:,1:28,5);
% avgShortNow = shortageCost(:,1:28,1,5);
figure;
subplot(3,1,1)
boxplot(avgMarNow, 'Labels',cellstr(string(s_P_abs(1:28))))
title('Distribution of Average MAR by mean P')
xlabel('Mean P [mm/m]')
ylabel('Mean MAR [MCM/y]')
subplot(3,1,2)
boxplot(avgStdNow ./ avgMarNow, 'Labels',cellstr(string(s_P_abs(1:28))))
title('Distribution of Average Coefficient of Variation by mean P')
xlabel('Mean P [mm/m]')
ylabel('Mean Std of Runoff [MCM/y]')
% subplot(3,1,3)
% boxplot(avgShortNow/1E6, 'Labels',cellstr(string(s_P_abs(1:28))))
% title('Distribution of Shortage Cost by mean P')
% xlabel('Mean P [mm/m]')
% ylabel('Shortage Cost M$')
% suptitle('Time Step 5: 2070-2090')

% Boxplots by T
avgMarNow = avgMAR(44:44+51,:,5)';
avgStdNow = avgStd(44:44+51,:,5)';
% avgShortNow = shortageCost(44:44+51,:,1,5)';
figure;
subplot(3,1,1)
boxplot(avgMarNow, 'Labels',cellstr(string(s_T_abs(44:44+51))))
title('Distribution of Average MAR by mean T')
xlabel('Mean T [degrees C]')
ylabel('Mean MAR [MCM/y]')
subplot(3,1,2)
boxplot(avgStdNow ./ avgMarNow, 'Labels',cellstr(string(s_T_abs(44:44+51))))
title('Distribution of Average Coefficient of Variation by mean T')
xlabel('Mean T [degrees C]')
ylabel('Mean Std of Runoff [MCM/y]')
% subplot(3,1,3)
% boxplot(avgShortNow/1E6, 'Labels',cellstr(string(s_T_abs(44:44+51))))
% title('Distribution of Shortage Cost by mean T')
% xlabel('Mean T [degrees C]')
% ylabel('Shortage Cost M$')
% suptitle('Time Step 5: 2070-2090')

% Add yield to these
end

%% Expansion policy for flexible dam

vldTInd = cell(1,5);
vldPInd = cell(1,5);
vldTInd{2} = 7:58;
vldPInd{2} = 9:16;
vldTInd{3} = 16:66;
vldPInd{3} = 4:22;
vldTInd{4} = 29:79;
vldPInd{4} = 3:23;
vldTInd{5} = 44:95;
vldPInd{5} = 1:28;

% Where there are single NAN gaps due to pruning, fill in with neighbor
XNow = cell(1,5);
XNow{2} = X(:,:,3,2);
XNow{2}(:,11) = XNow{2}(:,10);
XNow{2}(:,14) = XNow{2}(:,15);

XNow{3} = X(:,:,3,3);
XNow{3}(:,5) = XNow{3}(:,4);
XNow{3}(:,6) = XNow{3}(:,7);
XNow{3}(:,17) = XNow{3}(:,18);
XNow{3}(:,20) = XNow{3}(:,21);

XNow{4} = X(:,:,3,4);
XNow{4}(:,20) = XNow{4}(:,21);

XNow{5} = X(:,:,3,5);
% XNow{5}(52:54,:) = 0;

figure;
colormap( [0 0 0; .9 .9 .9])

for i = 2:5
    ax = subplot(2,2,i-1);
    imagesc([s_P_abs(vldPInd{i}(1)) s_P_abs(vldPInd{i}(end))],[s_T_abs(vldTInd{i}(1)) s_T_abs(vldTInd{i}(end))],  ...
        XNow{i}(vldTInd{i},vldPInd{i}))
    ax.YDir = 'normal';
%     xticklabels(cellstr(string(s_P_abs(vldPInd{i}))))
%     yticklabels(cellstr(string(s_T_abs(vldTInd{i}))))
    xlabel('Mean P [mm/m]')
    ylabel('Mean T [degrees C]')
    xlim([s_P_abs(1)-1 s_P_abs(28)+1])
    ylim([s_T_abs(1) s_T_abs(95)+.2])
    if i == 2
        hold on
        patch([100 100 101 101], [100 101 100 101], [0 0 0])
        patch([100 100 101 101], [100 101 100 101], [.9 .9 .9])
        legend('Do not expand','Expand')
        legend('location', 'NE')
    end
    title(decade{i})
end
suptitle('Expansion Policy for Flexible Dam')


%% Visualize value matrix for action 3 vs 4 in final time series

addpath(genpath('/Users/sarahfletcher/Documents/MATLAB/cbrewer'))


vldTInd{5} = 44:95;
vldPInd{5} = 1:28;

Vno = V(vldTInd{5},vldPInd{5},3,5);
Vexp = V(vldTInd{5},vldPInd{5},4,5);
Xexp = X(vldTInd{5},vldPInd{5},3,5);

figure;
ax = subplot(2,2,1);
imagesc([s_P_abs(vldPInd{5}(1)) s_P_abs(vldPInd{5}(end))],[s_T_abs(vldTInd{5}(1)) s_T_abs(vldTInd{5}(end))],  Vno )
colorbar
title('No expand V')
ax.YDir = 'normal';
ax = subplot(2,2,2);
imagesc([s_P_abs(vldPInd{5}(1)) s_P_abs(vldPInd{5}(end))],[s_T_abs(vldTInd{5}(1)) s_T_abs(vldTInd{5}(end))],  Vexp )
colorbar
title('Expand V')
ax.YDir = 'normal';
ax = subplot(2,2,3);
addpath(genpath('/Users/sarahfletcher/Documents/MATLAB/cbrewer'))
clrmp = cbrewer('div','RdYlBu',45);
colormap(clrmp)
imagesc([s_P_abs(vldPInd{5}(1)) s_P_abs(vldPInd{5}(end))],[s_T_abs(vldTInd{5}(1)) s_T_abs(vldTInd{5}(end))],  Vexp -Vno )
colorbar
title('Expand V - No expand V')
ax.YDir = 'normal';
ax = subplot(2,2,4);
imagesc([s_P_abs(vldPInd{5}(1)) s_P_abs(vldPInd{5}(end))],[s_T_abs(vldTInd{5}(1)) s_T_abs(vldTInd{5}(end))], Xexp)
colorbar
title('X')
ax.YDir = 'normal';




