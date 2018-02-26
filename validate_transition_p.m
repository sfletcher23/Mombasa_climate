%% Make T_clim and cum_T_clim


%% setup and make Transition mats
N = 9;
P_min = 50;
P_max = 84.5;
T_min = 22;
T_max = 30.8;
P_delta = 1.5; %mm/m
T_delta = 0.15; % deg C
s_P = P_min: P_delta: P_max;
s_T = T_min: T_delta: T_max;
T0 = 25.9;
P0 = 71;

[T_Precip, T_Temp] = samples2TClim(NUP(:,2:end,:), NUT(:,2:end,:), s_P, s_T, N);


%% Sample through T_Temp to get unconditional distribution
if false
samp = 10000;
p = rand(samp,N);
state_ind_P = zeros(samp,N);
state_ind_T = zeros(samp,N);
state_ind_P(:,1) = find(P0==s_P);
state_ind_T(:,1) = find(T0==s_T);
for i = 1:samp
    for t = 1:N-1
%         s_T(state_ind_T(i,t))
% %         state_ind_T(i,t)
%         p(i,t)
%          [s_T' T_Temp(:,state_ind_T(i,t),t)]
         
        state_ind_T(i,t+1) = find(p(i,t) < cumsum(T_Temp(:,state_ind_T(i,t),t)),1);
        state_ind_P(i,t+1) = find(p(i,t) < cumsum(T_Precip(:,state_ind_P(i,t),t)),1);
    end
end
T_over_time = s_T(state_ind_T);
P_over_time = s_P(state_ind_P);

T_over_time = s_T(state_ind_T);
P_over_time = s_P(state_ind_P);

figure;
s = randi(samp,10);
plot(T_over_time(s,:)')
ylabel('Decade')
xlabel('T')
title('10 samples of T samples over time')

mean_T_time = mean(T_over_time,1)

finalT = s_T(state_ind_T(:,end));
finalP = s_P(state_ind_P(:,end));

figure;
subplot(2,1,1)
hist(finalT, s_T)
title('Unconditional T in 2100: Tranistion probabilities')
xlim([27.4 29.2])
subplot(2,1,2)
hist(finalP,s_P)  
title('Unconditional P in 2100: Tranistion probabilities')
xlim([50 80])


T_p05 = prctile(T_over_time,.5);
T_p995 = prctile(T_over_time,99.5);
P_p05 = prctile(P_over_time,.5);
P_p995 = prctile(P_over_time,99.5);


[clrmp1]=cbrewer('seq', 'Reds', N);
[clrmp2]=cbrewer('seq', 'Blues', N);
figure;

t =1;
x = t:N;
X=[x,fliplr(x)];

subplot(2,1,1)
Y=[T_p05(t,t:end),fliplr(T_p995(t,t:end))];
hold on
fill(X,Y,clrmp1(t,:), 'LineWidth', 1.5);
xlabel('decade')
ylabel('degrees C')
title('99% CI for T time series')

subplot(2,1,2)
Y=[P_p05(t,t:end),fliplr(P_p995(t,t:end))];
hold on
fill(X,Y,clrmp2(t,:), 'LineWidth', 1.5); 
xlabel('decade')
ylabel('mm/m')
title('99% CI for P time series')


% figure;
% hist(NUT(:,2,find(s_T == 25)))
% title('T in decade 4 given 24.8 degrees in decade 3')

end
%% Updating over time
numSamp = 10000;
p = randi(numSamp,N-1);
state_ind_P = zeros(1,N);
state_ind_T = zeros(1,N);
state_ind_P(1) = find(P0==s_P);
state_ind_T(1) = find(T0==s_T);
T_over_time = cell(1,N);
P_over_time = cell(1,N);

for t = 1:N-1
    % Sample forward distribution given current state
    T_current = s_T(state_ind_T(t));
    P_current = s_P(state_ind_P(t));
    [T_over_time{t}, P_over_time{t}] = T2forwardSim(T_Temp, T_Precip, s_P, s_T, N, t, T_current, P_current, numSamp);

    % Sample next time period
    state_ind_T(t+1) = find(T_over_time{t}(p(t),2)==s_T);
    state_ind_P(t+1) = find(P_over_time{t}(p(t),2)==s_P);
end


% Lambdas
meanLamT_ts = zeros(21,N);
meanLamP_ts = zeros(21,N);
for t = 2:N
    meanLamT_ts(:,t) = mean(lamT(:,:,t,state_ind_T(t)) , 2);
    meanLamP_ts(:,t) = mean(lamP(:,:,t,state_ind_P(t)) , 2);
end
meanLamP_ts = meanLamP_ts(:,2:end);
meanLamT_ts = meanLamT_ts(:,2:end);
[clrmp3]=cbrewer('seq', 'Greys', 100);
[clrmp4]=cbrewer('seq', 'Reds', 100);

if true
figure
for t =1:N
    x = t:N;
    X=[x,fliplr(x)];
    T_p05 = prctile(T_over_time{t},.5);
    T_p995 = prctile(T_over_time{t},99.5);
    P_p05 = prctile(P_over_time{t},.5);
    P_p995 = prctile(P_over_time{t},99.5);
    
    subplot(2,2,1)
    Y=[T_p05,fliplr(T_p995)];
    hold on
    fill(X,Y,clrmp1(t,:), 'LineWidth', 1.5);
    scatter(t,s_T(state_ind_T(t)), 'k') 
    xlabel('decade')
    ylabel('degrees C')
    title('99% CI for T time series')
    
    subplot(2,2,2)
    colormap(clrmp3)
    imagesc(meanLamT_ts)
    colorbar
    xlabel('decade')
    ylabel('model')
    
    subplot(2,2,3)
    Y=[P_p05,fliplr(P_p995)];
    hold on
    fill(X,Y,clrmp1(t,:), 'LineWidth', 1.5);
    scatter(t,s_P(state_ind_P(t)), 'k') 
    xlabel('decade')
    ylabel('mm/m')
    title('99% CI for P time series')
    
    subplot(2,2,4)
    colormap(clrmp3)
    imagesc(meanLamP_ts)
    colorbar
    xlabel('decade')
    ylabel('model')
    
end
end

if false
figure
for t =1:N
    x = t:N;
    X=[x,fliplr(x)];
    T_p05 = prctile(T_over_time{t},.5);
    T_p995 = prctile(T_over_time{t},99.5);
    P_p05 = prctile(P_over_time{t},.5);
    P_p995 = prctile(P_over_time{t},99.5);
    
    subplot(1,2,1)
    Y=[T_p05,fliplr(T_p995)];
    hold on
    fill(X,Y,clrmp1(t,:), 'LineWidth', 1.5);
    scatter(t,s_T(state_ind_T(t)), 'k') 
    xlabel('decade')
    ylabel('degrees C')
    title('99% CI for T time series')
    
    subplot(1,2,2)
    colormap(clrmp3)
    imagesc(meanLamT_ts)
    colorbar
    xlabel('decade')
    ylabel('model')
   
end
end

%% Lambda plot

if false
% Lambda is models X samples X decades X states

[clrmp3]=cbrewer('seq', 'Blues', 100);

% Mean lamda for each model 
meanLamT_ts = zeros(21,N);
for t = 2:N
    meanLamT_ts(:,t) = mean(lamT(:,:,end,state_ind_T(t)) , 2);
end
meanLamT_ts = meanLamT_ts(:,2:end);
figure;
colormap(clrmp3)
imagesc(meanLamT_ts)
colorbar

% Normalized mean lambda for each model
norm_fct = sum(meanLamT_ts,1);
meanLamT_norm = meanLamT_ts ./ repmat(norm_fct,21,1);
figure;
colormap(clrmp3)
imagesc(meanLamT_norm)
colorbar
end

