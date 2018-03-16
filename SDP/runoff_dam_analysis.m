%% Runoff and dam sizing analysis


%% Analyze runoff 
if true

avgMAR = cellfun(@(x) mean(mean(x,2)), runoff);
avgStd = cellfun(@(x) mean(std(x,0,2)), runoff);


% Boxplots by P
avgMarNow = avgMAR(:,1:28,1);
avgStdNow = avgStd(:,1:28,1);
figure;
subplot(2,2,1)
boxplot(avgMarNow, 'Labels',cellstr(string(s_P_abs(1:28))))
title('Distribution of Average MAR by mean P')
xlabel('Mean P [mm/m]')
ylabel('Mean MAR [MCM/y]')
subplot(2,2,3)
boxplot(avgStdNow , 'Labels',cellstr(string(s_P_abs(1:28))))
title('Distribution of Avg Std Dev by mean P')
xlabel('Mean P [mm/m]')
ylabel('Mean Std of Runoff [MCM/y]')


% Boxplots by T
avgMarNow = avgMAR(44:44+51,:,1)';
avgStdNow = avgStd(44:44+51,:,1)';
subplot(2,2,2)
boxplot(avgMarNow, 'Labels',cellstr(string(s_T_abs(44:44+51))))
title('Distribution of Avg MAR by mean T')
xlabel('Mean T [degrees C]')
ylabel('Mean MAR [MCM/y]')
subplot(2,2,4)
boxplot(avgStdNow , 'Labels',cellstr(string(s_T_abs(44:44+51))))
title('Distribution of Avg Std Dev by mean T')
xlabel('Mean T [degrees C]')
ylabel('Mean Std of Runoff [MCM/y]')



%% Plot T and P time series

if true
    
    figure;
    for i = 3:12
        subplot(5,2,i-2)
        plot(P_ts{i,5}')
        hold on
        meanP = mean(mean(P_ts{i,5}',2));
        line([0 240], [meanP meanP], 'Color', 'k')
        ylim([0 400])
        xlim([0 240])
        title(strcat('P state: ', num2str(s_P_abs(i)), ' Mean P:', num2str(meanP)))
    end

end

%% Analyze yield and shortage costs

storage = [105 130];

unmet_ag = cell(M_T_abs, M_P_abs, length(storage), N);
unmet_dom = cell(M_T_abs, M_P_abs, length(storage), N);
yield = cell(M_T_abs, M_P_abs, length(storage), N);

for t = 1
    index_s_p_thisPeriod = index_s_p_time{t}; 
    for index_s_p = index_s_p_thisPeriod

        index_s_t_thisPeriod = index_s_t_time{t}; 
        for index_s_t= index_s_t_thisPeriod

            for s = 1:length(storage)

                [yield_mdl, K, dmd, unmet_dom_mdl, unmet_ag_mdl]  = ...
                    runoff2yield(runoff{index_s_t,index_s_p,t}, T_ts{index_s_t,t}, P_ts{index_s_p,t}, storage(s), runParam, climParam);
                yield{index_s_t, index_s_p, s, t} = yield_mdl;
                unmet_ag{index_s_t, index_s_p, s, t} = unmet_ag_mdl;
                unmet_dom{index_s_t, index_s_p, s, t} = unmet_dom_mdl;
%                 sum(unmet_dom_mdl,2)
%                 mean(sum(unmet_dom_mdl,2))
%                 storage(s)
            end
        end
    end
end




    
%% Sample single time series
if false

i = randi(125);
index_s_t = 44;
index_s_p = 1;
t = 5;

storage = [80 120 160];
figure;
for s = 1:length(storage)
    [yield_mdl, K, dmd, unmet_dom_mdl, unmet_ag_mdl]  = ...
        runoff2yield(runoff{index_s_t,index_s_p,t}(i:i+1,:), T_ts{index_s_t,t}(i:i+1,:), P_ts{index_s_p,t}(i:i+1,:), storage(s), runParam, climParam);
    subplot(length(storage), 1, s)
    plot(yield_mdl(1,:)')
    hold on
    plot(K(1,:)')
    plot(unmet_dom_mdl(1,:)')
    title(strcat('Unmet demand:', {' '}, num2str(sum(unmet_dom_mdl(1,:))), {' '}, 'Storage: ', {' '}, num2str(storage(s))))
end
legend('yield', 'storage', 'unmet')

end



%% Heatmaps

avgMAR = cellfun(@(x) mean(mean(x,2)), runoff);
avgP = cellfun(@(x) mean(mean(x,2)), P_ts);
avgT = cellfun(@(x) mean(mean(x,2)), T_ts);
avgUnmet = cellfun(@(x) mean(sum(x,2)), unmet_dom(:,:,1,:));

% Heatmap 
figure;
clrmp = cbrewer('div','RdYlBu',40);
colormap(clrmp)

subplot(2,2,1)
avgMarNow = avgMAR(44:44+51,1:28,1);
imagesc([s_P_abs(1) s_P_abs(28)],[s_T_abs(44) s_T_abs(44+51)], flipud(avgMarNow))
colorbar
xticks(s_P_abs(1): s_P_abs(28))
yticks(s_T_abs(44):.5: s_T_abs(44+51))
xticklabels(s_P_abs(1:28))
yticklabels(fliplr(s_T_abs(44):.5:s_T_abs(44+51)))
xlabel('P [mm/m]')
ylabel('T [degrees C]')
title('Mean Runoff')

subplot(2,2,2)
colormap(gca, flipud(clrmp))
avgUnmetNow = avgUnmet(44:44+51,1:28,1,1);
imagesc([s_P_abs(1) s_P_abs(28)],[s_T_abs(44) s_T_abs(44+51)], flipud(avgUnmetNow))
colorbar
xticks(s_P_abs(1): s_P_abs(28))
yticks(s_T_abs(44):.5: s_T_abs(44+51))
xticklabels(s_P_abs(1:28))
yticklabels(fliplr(s_T_abs(44):.5:s_T_abs(44+51)))
xlabel('P [mm/m]')
ylabel('T [degrees C]')
title('Mean Unmet Demand')

subplot(2,2,3)
colormap(gca, clrmp)
avgPNow = avgP(1:28,1);
imagesc([s_P_abs(1) s_P_abs(28)],[s_T_abs(44) s_T_abs(44+51)],avgPNow')
colorbar
yticks([])
xticks(s_P_abs(1): s_P_abs(28))
xticklabels(s_P_abs(1:28))
xlabel('P [mm/m]')
title('Mean P')

subplot(2,2,4)
colormap(gca, flipud(clrmp))
avgTNow = avgT(44:44+51,1);
imagesc([s_P_abs(1) s_P_abs(28)],[s_T_abs(44) s_T_abs(44+51)],flipud(avgTNow))
colorbar
xticks([])
yticks(s_T_abs(44):.5: s_T_abs(44+51))
yticklabels(fliplr(s_T_abs(44):.5:s_T_abs(44+51)))
ylabel('T [degrees C]')
title('Mean T')


end

%% Plot shortage 

avgTotalUnmet = cellfun(@(x) mean(sum(x,2)), unmet_dom(:,:,end,1));
p90TotalUnmet = cellfun(@(x) prctile(sum(x,2),90), unmet_dom(:,:,end,1));
p10TotalUnmet = cellfun(@(x) prctile(sum(x,2),10), unmet_dom(:,:,end,1));

% Boxplots by P
avgUnmetNow = avgTotalUnmet(:,1:28);
p10UnmetNow = p10TotalUnmet(:,1:28);
p90UnmetNow = p90TotalUnmet(:,1:28);
figure;
subplot(3,2,1)
boxplot(p10UnmetNow, 'Labels',cellstr(string(s_P_abs(1:28))))
title('p10 UnmetDemand by mean P')
xlabel('Mean P [mm/m]')
ylabel('Unmet demand [MCM]')
subplot(3,2,3)
boxplot(avgUnmetNow, 'Labels',cellstr(string(s_P_abs(1:28))))
title('Avg UnmetDemand by mean P')
xlabel('Mean P [mm/m]')
ylabel('Unmet demand [MCM]')
subplot(3,2,5)
boxplot(p90UnmetNow, 'Labels',cellstr(string(s_P_abs(1:28))))
title('p90 Unmet Demand by mean P')
xlabel('Mean P [mm/m]')
ylabel('Unmet demand [MCM]')


% Boxplots by T
avgUnmetNow = avgTotalUnmet(44:44+51,:)';
p10UnmetNow = p10TotalUnmet(44:44+51,:)';
p90UnmetNow = p90TotalUnmet(44:44+51,:)';
subplot(3,2,2)
boxplot(p10UnmetNow, 'Labels',cellstr(string(s_T_abs(44:44+51))))
title('p10 UnmetDemand by mean T')
xlabel('Mean T [degrees C]')
ylabel('Unmet demand [MCM]')
subplot(3,2,4)
boxplot(avgUnmetNow, 'Labels',cellstr(string(s_T_abs(44:44+51))))
title('Avg UnmetDemand by mean T')
xlabel('Mean T [degrees C]')
ylabel('Unmet demand [MCM]')
subplot(3,2,6)
boxplot(p90UnmetNow, 'Labels',cellstr(string(s_T_abs(44:44+51))))
title('p90 UnmetDemand by mean T')
xlabel('Mean T [degrees C]')
ylabel('Unmet demand [MCM]')

suptitle('Storage 80 MCM, BiDecade 5')


%% Plot yield and shortage costs 

avgYield = cellfun(@(x) mean(prctile(x,10,2)), yield(:,:,end,1));
avgUnmet = cellfun(@(x) mean(prctile(x,90,2)), unmet_dom(:,:,end,1));


% Boxplots by P
avgYieldNow = avgYield(:,1:28);
avgUnmetNow = avgUnmet(:,1:28);
figure;
subplot(2,2,1)
boxplot(avgYieldNow, 'Labels',cellstr(string(s_P_abs(1:28))))
title('Distribution of p10 Yield by mean P')
xlabel('Mean P [mm/m]')
ylabel('Yield [MCM/y]')
subplot(2,2,3)
boxplot(avgUnmetNow, 'Labels',cellstr(string(s_P_abs(1:28))))
title('Distribution of p90 Unmet Demand by mean P')
xlabel('Mean P [mm/m]')
ylabel('Annual demand [MCM/y]')


% Boxplots by T
avgYieldNow = avgYield(44:44+51,:)';
avgUnmetNow = avgUnmet(44:44+51,:)';
subplot(2,2,2)
boxplot(avgYieldNow, 'Labels',cellstr(string(s_T_abs(44:44+51))))
title('Distribution of p10 Yield by mean T')
xlabel('Mean T [degrees C]')
ylabel('Yield [MCM/y]')
subplot(2,2,4)
boxplot(avgUnmetNow, 'Labels',cellstr(string(s_T_abs(44:44+51))))
title('Distribution of p90 Unmet Demand by mean T')
xlabel('Mean T [degrees C]')
ylabel('Annual demand [MCM/y]')

suptitle('Storage 140 MCM, BiDecade 5')

%% Plot umnmet demand by dam size

avgTotalUnmet = squeeze (cellfun(@(x) mean(sum(x,2)), unmet_dom(:,:,:,1)) );

indP = [1:8:32];
indT = 44+51;
figure;
for i = 1:4
    subplot(2,2,i)
    boxplot(squeeze(avgTotalUnmet(:,indP(i),:)))
    xticklabels(cellstr(num2str(storage')))
    ylim([0 1000])
    xlabel('Storage [MCM]')
    ylabel('Mean unmet demand [MCM]')
    title(strcat('P state: ', num2str(s_P_abs(indP(i)))))
end
suptitle('Unmet demand by dam storage (variation across T)')

indP = 16;
indT = 54:10:90;
figure;
for i = 1:4
    subplot(2,2,i)
    boxplot(squeeze(avgTotalUnmet(indT(i),:,:)))
    xticklabels(cellstr(num2str(storage')))
    ylim([0 1000])
    xlabel('Storage [MCM]')
    ylabel('Mean unmet demand [MCM]')
    title(strcat('T state: ', num2str(s_T_abs(indT(i)))))
end
suptitle('Unmet demand by dam storage (variation across P)')

%% Heatmaps Unmet demand by storage

avgUnmet = cellfun(@(x) mean(sum(x,2)), unmet_dom(:,:,:,:));


% Heatmap 
f = figure;
addpath(genpath('/Users/sarahfletcher/Documents/MATLAB/cbrewer'))
clrmp = cbrewer('div','RdYlBu',45);
colormap(clrmp)
font_size = 13;

for i = 1:2
    h(i) = subplot(3,2,i);
    colormap(gca, flipud(clrmp))
    avgUnmetNow = avgUnmet(44:44+51,1:28,i,1);
    imagesc([s_P_abs(1) s_P_abs(28)],[s_T_abs(44) s_T_abs(44+51)], flipud(avgUnmetNow))
    c = colorbar;
    caxis([0 300])
    c.Label.String = 'MCM over 20 yrs';
    xticks(s_P_abs(1):2: s_P_abs(28))
    yticks(28:.5:30.5)
    xticklabels(s_P_abs(1:2:28))
    yticklabels(sprintf('%.1f\n', 30.5:-.5:28))
    xlabel('P [mm/m]')
    ylabel('T [degrees C]')
    title(strcat('Storage: ', {' '}, num2str(storage(i)), {' '}, 'MCM'))
end
suptitle('Unmet Demand by Reservoir Size')
set(findall(f,'-property','FontSize'),'FontSize',font_size)


