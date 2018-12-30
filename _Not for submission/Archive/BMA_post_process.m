%% BMA post processing


% NU and MU matrices are numSamples (10000) x numDecades (10) x numScenarios (3)

%% Trace plots

decades = {'2000-2009', '2010-2019', '2020-2029', '2030-2039', '2040-2049', ...
    '2050-2059', '2060-2069', '2070-2079', '2080-2089', '2090-2099'};
N = 10;
numScenarios = 59;
saveOn = false;

% Future precip
if false
for j = 1:numScenarios % make a new figure for each scenario
    figure;
    for i = 1:10 % make a new subplot for each decade
        subplot(2,5,i)
        plot(1:1000,NUP(:,i,j))
        title(decades{i})
        ylim([0, 100])
        if i == 1 || i == 6
            ylabel('P (mm)')
        end
        if i >= 6
            xlabel('iteration')
        end
    end
    ttl = strcat('Future P Trace Plot: Scenario:', ' ', num2str(j));
    suptitle(ttl)
    if saveOn
        savefig(gcf, strcat('BMA_analysis_', date, '/', ttl))
    end
end
end

% Future temp
if true
for j = 1:numScenarios % make a new figure for each scenario
    figure;
    for i = 1:10 % make a new subplot for each decade
        subplot(2,5,i)
        plot(1:1000,NUT(:,i,j))
        title(decades{i})
        ylim([23, 31])
        if i == 1 || i == 6
            ylabel('T (degrees C)')
        end
        if i >= 6
            xlabel('iteration')
        end
    end
    ttl = strcat('Future T Trace Plot: Scenario:', ' ', num2str(j));
    suptitle(ttl)
    if saveOn
        savefig(strcat('BMA_analysis_', date, '/', ttl))
    end
end
end

% Hist precip
if false
for j = 1:numScenarios % make a new figure for each scenario
    figure;
    for i = 1:10 % make a new subplot for each decade
        subplot(2,5,i)
        plot(1:1000,MUP(:,i,j))
        title(decades{i})
        ylim([0, 100])
        if i == 1 || i == 6
            ylabel('P (mm)')
        end
        if i >= 6
            xlabel('iteration')
        end
    end
    ttl = strcat('Hist P Trace Plot: Scenario:', ' ', num2str(j));
    suptitle(ttl)
    if saveOn
        savefig(strcat('BMA_analysis_', date, '/', ttl))
    end
end
end

% Hist temp
if false
for j = 1:numScenarios % make a new figure for each scenario
    figure;
    for i = 1:10 % make a new subplot for each decade
        subplot(2,5,i)
        plot(1:1000,MUT(:,i,j))
        title(decades{i})
        ylim([20, 31])
        if i == 1 || i == 6
            ylabel('T (degrees C)')
        end
        if i >= 6
            xlabel('iteration')
        end
    end
    ttl = strcat('Hist T Trace Plot: Scenario:', ' ', num2str(j));
    suptitle(ttl)
    if saveOn
        savefig(strcat('BMA_analysis_', date, '/', ttl))
    end
end
end

%% Trace plot zoom
% Comparing last decade with other decade at higher resolution to assess
% cycle problem

% Future Precip
sampleIndex = 900:1000; % Using a subset of the samples for better visualization
j = 1; % Scenario 
i = 10; % Decade 1
k = 9; % Deacde 2
figure;
subplot(2,1,1)
plot(sampleIndex, NUP(sampleIndex,i,j),'-o')
title(decades{i})
subplot(2,1,2)
plot(sampleIndex, NUP(sampleIndex,k,j),'-o')
title(decades{k})
ttl = strcat('Zoomed trace plot fut P scenario:', num2str(j));
suptitle(ttl)
if saveOn
    savefig(strcat('BMA_analysis_', date, '/', ttl))
end
% I played around with looking at different decades and scenarios - can
% really see the 21-step cycle. Although it doesn't draw the same exaclty values each time -
% suggests it is cycling through the models in a purely deterministic
% fashion but sampling with variation from each. Can we look at the marginal samples
% from another parameter, maybe lambda?, to confirm this? 

% Future Temp
sampleIndex = 900:1000; % Using a subset of the samples for better visualization
j = 1; % Scenario 5
figure;
subplot(2,1,1)
plot(sampleIndex, NUT(sampleIndex,i,j),'-o')
title(decades{i})
subplot(2,1,2)
plot(sampleIndex, NUT(sampleIndex,k,j),'-o')
title(decades{k})
ttl = strcat('Zoomed trace plot fut T scenario 5:', num2str(j));
suptitle(ttl)
if saveOn
    savefig(strcat('BMA_analysis_', date, '/', ttl))
end
% Can clearly see the difference between cycle and non-cycle


%% Autocorrelation

step = 1; % Increase the step to discard more of the samples
sampleIndex = 1:step:1000; 

% Future precip
if false
for j = 1:numScenarios % make a new figure for each scenario
    figure;
    for i = 1:10 % make a new subplot for each decade
        subplot(2,5,i)
        autocorr(NUP(sampleIndex,i,j),44)
        title(decades{i})
    end
    ttl = strcat('Future P Autocorrelation: Scenario:', ' ', num2str(j));
    suptitle(ttl)
    if saveOn
        savefig(strcat('BMA_analysis_', date, '/', ttl))
    end
end
end

% Future temp
if true
for j = 1:numScenarios % make a new figure for each scenario
    figure;
    for i = 1:10 % make a new subplot for each decade
        subplot(2,5,i)
        autocorr(NUT(sampleIndex,i,j), 44)
        title(decades{i})
    end
    ttl = strcat('Future T Autocorrelation: Scenario:', ' ', num2str(j));
    suptitle(ttl)
    if saveOn
        savefig(strcat('BMA_analysis_', date, '/', ttl))
    end
end
end

%% Plot future distributions

% Temp
figure;
deltaT = 1:10:51;
c = 1;
for i = 1:10:51
    subplot(6, 1, c)
    boxplot(NUT(:,:,i))
    ylim([24 31]) 
    xticklabels(decades)
    txt = 'Total Delta T: %.2f';
    %title(sprintf(txt, deltaT(i)));
    ylabel('Degrees C')
    c = c+1;
end

figure;
deltaP = 1:5:21
c = 1;
for i = 1:5:21
    subplot(5, 1, c)
    boxplot(NUP(:,:,i))
    ylim([45 90]) 
    xticklabels(decades)
    txt = 'Total Delta P: %.2f %%';
    %title(sprintf(txt, s(i)));
    ylabel('mm/m')
    c = c+1;
end

%% Histograms to assess discretization

figure; 
subplot(1,4,1)
hist(NUT(:,9,3), 24:0.05:30)
title('0.05 C')
subplot(1,4,2)
hist(NUT(:,9,3), 24:0.1:30)
title('0.1 C')
subplot(1,4,3)
hist(NUT(:,9,3), 24:0.15:30)
title('0.15 C')
subplot(1,4,4)
hist(NUT(:,9,3), 24:0.2:30)
title('0.2 C')
suptitle('T discretization ')

figure; 
subplot(1,4,1)
hist(NUP(:,9,1), 50:0.5:85)
title('0.5 mm/d')
subplot(1,4,2)
hist(NUP(:,9,1), 50:1:85)
title('1 mm/d')
subplot(1,4,3)
hist(NUP(:,9,1), 50:1.5:85)
title('1.5 mm/d')
subplot(1,4,4)
hist(NUP(:,9,1), 50:2:85)
title('2 mm/d')
suptitle('P discretization ')

%% Make T_clim and cum_T_clim

N = 10;
P_min = 50;
P_max = 84.5;
T_min = 23.6;
T_max = 30.8;
P_delta = 1.5; %mm/d
T_delta = 0.15; % deg C
s_P = P_min: P_delta: P_max;
s_T = T_min: T_delta: T_max;

[T_Precip, T_Temp] = samples2TClim(NUP, NUT, s_P, s_T, N);

% Sample through T_Temp to get unconditional distribution
samp = 100000;
T0 = s_T(1);
P0 = s_P(12);
p = rand(samp,N);
state_ind_P = zeros(samp,N);
state_ind_T = zeros(samp,N);
state_ind_P(:,1) = find(P0==s_P);
state_ind_T(:,1) = find(T0==s_T);
for i = 1:samp
    for t = 1:N-1
        state_ind_T(i,t+1) = find(p(i,t) < cumsum(T_Temp(:,state_ind_T(t),t)),1);
        state_ind_P(i,t+1) = find(p(i,t) < cumsum(T_Precip(:,state_ind_P(t),t)),1);
    end
end
T_over_time = s_T(state_ind_T);
P_over_time = s_P(state_ind_P);

figure;
s = randi(samp,10);
plot(T_over_time(s,:)')
mean_T_time = mean(s_T(state_ind_T),1);

finalT = s_T(state_ind_T(:,end));
finalP = s_P(state_ind_P(:,end));



figure;
subplot(2,1,1)
hist(finalT, s_T)
subplot(2,1,2)
hist(finalP,s_P)    


% 