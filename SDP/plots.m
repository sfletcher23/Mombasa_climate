%% Analyze runoff and shortage

avgMAR = cellfun(@(x) mean(mean(x,2)), runoff);
avgStd = cellfun(@(x) mean(std(x,0,2)), runoff);
short = unmet_ag + unmet_dom;

%%

% Boxplots by P
avgMarNow = avgMAR(:,1:28,5);
avgStdNow = avgStd(:,1:28,5);
avgShortNow = short(:,1:28,1,5);
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
subplot(3,1,3)
boxplot(avgShorNow, 'Labels',cellstr(string(s_P_abs(1:28))))
title('Distribution of Shortage by mean P')
xlabel('Mean P [mm/m]')
ylabel('Mean Std of Runoff [MCM/y]')
suptitle('Time Step 5: 2070-2090')

% Boxplots by T
avgMarNow = avgMAR(44:44+51,:,5)';
avgStdNow = avgStd(44:44+51,:,5)';
figure;
subplot(2,1,1)
boxplot(avgMarNow, 'Labels',cellstr(string(s_T_abs(44:44+51))))
title('Distribution of Average MAR by mean T')
xlabel('Mean T [degrees C]')
ylabel('Mean MAR [MCM/y]')
subplot(2,1,2)
boxplot(avgStdNow ./ avgMarNow, 'Labels',cellstr(string(s_T_abs(44:44+51))))
title('Distribution of Average Coefficient of Variation by mean T')
xlabel('Mean T [degrees C]')
ylabel('Mean Std of Runoff [MCM/y]')
suptitle('Time Step 5: 2070-2090')
