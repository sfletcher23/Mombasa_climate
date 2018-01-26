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
figure;
for i = 1:numDecades
    subplot(5, 2, i)
    mean_P = mean(mean(P_ts(:,:,i,4)));
    plot(1:numMonths, P_ts(:,1:10,i,4))
    ylim([0 400])
    title(strcat(decades{i}, {' '}, 'Mean P:', {' '}, num2str(mean_P, '%.1f')))
end

figure;
for i = 1:numDecades
    subplot(5, 2, i)
    mean_T = mean(mean(T_ts(:,:,i,4)));
    plot(1:numMonths, T_ts(:,1:10,i,4))
    ylim([22 30])
    title(strcat(decades{i}, {' '}, 'Mean T:', {' '}, num2str(mean_T, '%.1f')))
end

%% Use ClIRUN to calculate streamflow

    
    
    
    