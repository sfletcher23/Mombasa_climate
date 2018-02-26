load('BMA_data_1-26-18')

% P_ts and T_ts are 4d matrices with monthly 10-year time series samples of
% P and T respectively. They are numMonths (120) x numSamples (100) x
% numDecades (10) x numScenario (9)

[numMonths, numSamples, numDecades, numScenarios] = size(P_ts);

% Scenarios
% for T: 
% Scenario 1: no temp change
% Scenario 2: temp change of 2.25 by end of century
% Scenario 3: temp change of 4.5 by end of century
% 
% for P: 
% Scenario 1: decrease of P by 72% by end of century
% Scenario 2: no change
% Scenario 3: increase of P by 72% by end of century
% 
% Scenario order:
% P1, T1
% P1, T2
% P1, T3
% P2, T1  <-- this one is no change
% P2, T2
% P2, T3
% P3, T1
% P3, T2
% P3, T3

decades = {'2000-2009', '2010-2019', '2020-2029', '2030-2039', '2040-2049', ...
    '2050-2059', '2060-2069', '2070-2079', '2080-2089', '2090-2099'};

%% Plot sample time series, no change (scen4)

if true
    figure;
    for i = 1:numDecades
        subplot(5, 2, i)
        mean_P = mean(mean(P_ts(:,:,i,4)));
        std_P = std(mean(P_ts(:,:,i,4)));
        plot(1:numMonths, P_ts(:,1:10,i,4))
        ylim([0 400])
        title(strcat(decades{i}, {' '}, 'Mean P:', {' '}, num2str(mean_P, '%.1f') ...
            , {' '}, 'Std:', {' '}, num2str(std_P, '%.1f')))
    end

    figure;
    for i = 1:numDecades
        subplot(5, 2, i)
        mean_T = mean(mean(T_ts(:,:,i,4)));
        std_T = std(mean(T_ts(:,:,i,4)));
        plot(1:numMonths, T_ts(:,1:10,i,4))
        ylim([22 30])
        title(strcat(decades{i}, {' '}, 'Mean T:', {' '}, num2str(mean_T, '%.1f') ...
            , {' '}, 'Std:', {' '}, num2str(std_T, '%.1f')))
    end
end
%% Use ClIRUN to calculate streamflow
calibrationFile = '16_Jan_2018_14_39_45_ogdata_point_normparam_1';
calibrationFile = '29_Jan_2018_17_10_19_maybe_winner_3.mat';
load(strcat('CLIRUN/OutputData/data/',calibrationFile), 'X_results')

streamflow_mmpd = zeros(numSamples,numMonths);
for i = 1:numSamples
    T = T_ts(:,i,1,4)';
    P = P_ts(:,i,1,4)';
    watyear = 1;
    streamflow_mmpd(i,:) = Simulator(X_results, T, P, watyear); %mm/d
end

area = 2250 * 1E6; %m2
streamflow_cmpd = streamflow_mmpd /1E3 * area;
streamflow_mcmpy = cmpd2mcmpy(streamflow_cmpd);
figure;
plot(streamflow_mcmpy(1:10,:)')   
mar_mcmpy = mean(mean(streamflow_mcmpy));
%% Compare estimates for firm yield and yield from current no change scenario to historical

% Firm yield
storage = 120;
fyield = zeros(numSamples,1);
for i = 1:numSamples
    inflow = streamflow_mcmpy(i,:);
    fyield(i) = strmflw2frmyld(inflow, storage);
end
avg_fyield = mean(fyield)


load('Mwache_hydro_data', 'runoff_mwache', 'P0_date_lin', 'T0_date_lin')
hist_strflw = runoff_mwache /1E3 * area * 365 / 1E6; % convert from mm/d to mcm/y
hist_fyield = strmflw2frmyld(hist_strflw, storage)

% Yield simulations

[yield_mdl, K, dmd, unmet_dom_mdl, unmet_ag_mdl]  = inflow2yield(hist_strflw, T0_date_lin', P0_date_lin', storage);
figure; 
subplot(2,1,2)
bar([yield_mdl; unmet_dom_mdl; unmet_ag_mdl]', 'stacked');
hold on
plot(dmd, 'LineWidth', 1.5)
xlim([0 length(hist_strflw)])
legend('Yield', 'Unmet domestic', 'Unmet ag', 'Demand')
subplot(2,1,1)
hold on
plot(K+20, 'LineWidth', 1.5)
plot(hist_strflw, 'LineWidth', 1.5)
legend('storage', 'inflow')
title('historical inflow')
xlim([0 length(hist_strflw)])


i = randi(numSamples);
[yield_mdl, K, dmd, unmet_dom_mdl, unmet_ag_mdl]  = inflow2yield(streamflow_mcmpy(i,:), T_ts(:,i,1,4)', P_ts(:,i,1,4)', storage);
figure; 
subplot(2,1,2)
bar([yield_mdl; unmet_dom_mdl; unmet_ag_mdl]', 'stacked');
hold on
plot(dmd, 'LineWidth', 1.5)
xlim([0 length(streamflow_mcmpy(i,:))])
legend('Yield', 'Unmet domestic', 'Unmet ag', 'Demand')
subplot(2,1,1)
hold on
plot(K+20, 'LineWidth', 1.5)
plot(streamflow_mcmpy(i,:), 'LineWidth', 1.5)
legend('storage', 'inflow')
title('simulated inflow')
xlim([0 length(streamflow_mcmpy(i,:))])



