%% BMA samples to transition probabilities and forward simulation of T and P




%% Setup and make delta Transition mats

if true
load('BMA_results_deltap05T_p2P07-Feb-2018 20:18:49.mat')
end

N = 5;
P_min = -.3;
P_max = .3;
P_delta = .02; %mm/m
s_P = P_min:P_delta:P_max;
P0 = s_P(15);
P0_abs = 75;
M_P = length(s_P);

T_min = 0;
T_max = 1.5;
T_delta = 0.05; % deg C
s_T = T_min: T_delta: T_max;
T0 = s_T(1);
T0_abs = 26;
M_T = length(s_T);

decades = { '1990', '2010', '2030', '2050', '2070', '2090'};

[ T_Temp] = samples2TTemp( NUT, s_T, N);
[ T_Precip] = samples2TTemp( NUP, s_P, N);

numSamp = 10000;

%% Sample through T_Temp to get delta time series

if true

p = rand(numSamp,N+1);
state_ind_P = zeros(numSamp,N);
state_ind_T = zeros(numSamp,N+1);
state_ind_P(:,1) = find(P0==s_P);
T0_ind = randi(M_T,numSamp,1);

state_ind_T(:,1) = T0_ind;
for i = 1:numSamp
    for t = 1:N
        state_ind_T(i,t+1) = find(p(i,t) < cumsum(T_Temp(:,state_ind_T(i,t),t)),1);
        state_ind_P(i,t+1) = find(p(i,t) < cumsum(T_Precip(:,state_ind_P(i,t),t)),1);
    end
end
T_delta_over_time = s_T(state_ind_T);
T_p05 = prctile(T_delta_over_time,.5);
T_p995 = prctile(T_delta_over_time,99.5);
P_delta_over_time = s_P(state_ind_P);
P_p05 = prctile(P_delta_over_time,.5);
P_p995 = prctile(P_delta_over_time,99.5);

[clrmp1]=cbrewer('seq', 'Reds', N);
[clrmp2]=cbrewer('seq', 'Blues', N);

if false
figure;
t =1;
x = t:N+1;
X=[x,fliplr(x)];

subplot(2,1,1)
Y=[T_p05(t,t:end),fliplr(T_p995(t,t:end))];
hold on
fill(X,Y,clrmp1(t,:), 'LineWidth', 1.5);
xlabel('decade')
ylabel('delta degrees C')
title('99% CI for delta T time series')

subplot(2,1,2)
Y=[P_p05(t,t:end),fliplr(P_p995(t,t:end))];
hold on
fill(X,Y,clrmp1(t,:), 'LineWidth', 1.5);
xlabel('decade')
ylabel('delta degrees C')
title('99% CI for delta T time series')
end
end

%% Calculate absolute values from deltas

if true

% Sum delta to get absolutes
T_over_time = cumsum( T_delta_over_time,2) + T0_abs;

% Precip is percent change
P_over_time = zeros(numSamp, N+1);
for i = 1:numSamp
    for t = 1:6
        if t == 1 
            P_over_time(i,t) = P0_abs;
        else
            P_over_time(i,t) = P_over_time(i,t-1) * (1+P_delta_over_time(i,t));
        end
    end
end


% Plots
if false
    
figure;
s = randi(numSamp,10);
subplot(1,2,1)
plot(T_over_time(s,:)')
xlabel('Decade')
ylabel('T')
title('10 numSamples of T numSamples over time')
subplot(1,2,2)
plot(P_over_time(s,:)')
xlabel('Decade')
ylabel('P')
title('10 numSamples of P numSamples over time')

mean_T_time = mean(T_over_time,1);
mean_P_time = mean(P_over_time,1);

finalT = s_T(state_ind_T(:,end));
finalP = s_P(state_ind_P(:,end));

figure;
subplot(2,1,1)
hist(finalT, s_T, 100)
title('Unconditional T in 2100: Tranistion probabilities')
%xlim([27.4 29.2])
subplot(2,1,2)
hist(finalP,s_P)  
title('Unconditional P in 2100: Tranistion probabilities')
% xlim([50 80])

T_p05 = prctile(T_over_time,.5);
T_p995 = prctile(T_over_time,99.5);
P_p05 = prctile(P_over_time,.5);
P_p995 = prctile(P_over_time,99.5);


[clrmp1]=cbrewer('seq', 'Reds', N+1);
[clrmp2]=cbrewer('seq', 'Blues', N+1);
figure;

t =1;
x = t:N+1;
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
end


end

%% Use time series to caluclate absolute transition probabilities 

T_abs_max = max(s_T) * N;
s_T_abs = T0_abs:T_delta:T0_abs+ T_abs_max;
M_T_abs = length(s_T_abs);
T_bins = [s_T_abs-T_delta/2 s_T_abs(end)+T_delta/2];
T_Temp_abs = zeros(M_T_abs,M_T_abs,N);

P_abs_max = max(s_P) * N;
s_P_abs = 66:1:97;
M_P_abs = length(s_P_abs);
P_bins = [s_P_abs-P_delta/2 s_P_abs(end)+P_delta/2];
T_Precip_abs = zeros(M_P_abs,M_P_abs,N);

for i = 1:length(s_T_abs)
    for t = 1:N
        T_current = s_T_abs(i);
        indexNow = find(T_over_time(:,t) == T_current);
        relevant_deltas = T_over_time(indexNow,t+1) - T_over_time(indexNow,t);
        T_next = T_current + relevant_deltas;
        T_Temp_abs(:,i,t) = histcounts(T_next, T_bins, 'Normalization', 'Probability');     
    end
end

P_over_time_rounded = round(P_over_time);
for i = 1:length(s_P_abs)
    for t = 1:N
        P_current = s_P_abs(i);
        indexNow = find(P_over_time_rounded(:,t) == P_current);
        relevant_deltas = P_over_time_rounded(indexNow,t+1) - P_over_time_rounded(indexNow,t);
        P_next = P_current + relevant_deltas;
        T_Precip_abs(:,i,t) = histcounts(P_next, P_bins, 'Normalization', 'Probability');     
    end
end

%% Use abs trans ps to simulate time series
if false
[T_over_time_abs] = T2forwardSimTemp(T_Temp_abs, s_T_abs, N, 1, T0_abs, 10000, true);
figure;
plot(1:N+1,T_over_time_abs')

[P_over_time_abs] = T2forwardSimTemp(T_Precip_abs, s_P_abs, N, 1, s_P_abs(133), 10000, false);
figure;
plot(1:N+1,P_over_time_abs')

% 
% numSamp = 100000;
% p = rand(numSamp,N+1);
% state_ind_T = zeros(numSamp,N+1);
% state_ind_T(:,1) = find(T0_abs==s_T_abs);
% for i = 1:numSamp
%     for t = 1:N
%         state_ind_T(i,t+1) = find(p(i,t) < cumsum(T_Temp_abs(:,state_ind_T(i,t),t)),1);
%     end
% end
% T_over_time_abs = s_T(state_ind_T);
% T_p05 = prctile(T_delta_over_time_abs,.5);
% T_p995 = prctile(T_delta_over_time_abs,99.5);

end





%% Updating over time
if true


p = randi(numSamp,N-1);
state_ind_P = zeros(1,N);
state_ind_T = zeros(1,N);
state_ind_P(1) =  find(P0_abs==s_P_abs);
state_ind_T(1) = find(T0_abs==s_T_abs);
% T0_ind = randi(31,1,1);
% state_ind_T(1) = T0_ind;
T_over_time = cell(1,N);
P_over_time = cell(1,N);

for t = 1:N
    % Sample forward distribution given current state
    T_current = s_T_abs(state_ind_T(t));
    P_current = s_P_abs(state_ind_P(t));
    [T_over_time{t}] = T2forwardSimTemp(T_Temp_abs, s_T_abs, N, t, T_current, numSamp, false);
    [P_over_time{t}] = T2forwardSimTemp(T_Precip_abs, s_P_abs, N, t, P_current, numSamp, false);

    
    % Sample next time period
    state_ind_T(t+1) = find(T_over_time{t}(p(t),2)==s_T_abs);
    state_ind_P(t+1) = find(P_over_time{t}(p(t),2)==s_P_abs);
end


% Lambdas
if false
[clrmp3]=cbrewer('seq', 'Greys', 100);
[clrmp4]=cbrewer('seq', 'Reds', 100);
states = 1:5:30;
figure
for i = 1:length(states)
    meanLamT_ts = zeros(21,N);
    meanLamP_ts = zeros(21,N);
    for t = 1:N
        meanLamT_ts(:,t) = mean(lamT(:,:,t,states(i)) , 2);
        meanLamP_ts(:,t) = mean(lamP(:,:,t,states(i)), 2);
    end
    meanLamP_ts = meanLamP_ts(:,2:end);
    meanLamT_ts = meanLamT_ts(:,2:end);
    subplot(2,3,i)
    colormap(clrmp3)
    imagesc(meanLamT_ts, [0 100])
    colorbar
    xlabel('decade')
    ylabel('model')
    title(strcat('State: ', {' '}, num2str(s_T(i))))
end
end


if false
[clrmp3]=cbrewer('seq', 'Greys', 100);
[clrmp4]=cbrewer('seq', 'Reds', 100);
figure
for t =1:N
    x = t:N+1;
    X=[x,fliplr(x)];
    T_p05 = prctile(T_over_time{t},.01);
    T_p995 = prctile(T_over_time{t},99.9);
    P_p05 = prctile(P_over_time{t},.01);
    P_p995 = prctile(P_over_time{t},99.9);
    
    subplot(2,2,1)
    Y=[T_p05,fliplr(T_p995)];
    hold on
    fill(X,Y,clrmp1(t,:), 'LineWidth', 1.5);
    scatter(t,s_T_abs(state_ind_T(t)), 'k') 
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
    scatter(t,s_P_abs(state_ind_P(t)), 'k') 
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

if true

  
fig = figure;
[clrmp1]=cbrewer('seq', 'Reds', N);
[clrmp2]=cbrewer('seq', 'Blues', N);
set(fig,'Position', [680 558 1400 750])


for t =1:N
    x = t:N+1;
    X=[x,fliplr(x)];
    T_p05 = prctile(T_over_time{t},.01);
    T_p995 = prctile(T_over_time{t},99.9);
    P_p05 = prctile(P_over_time{t},.01);
    P_p995 = prctile(P_over_time{t},99.9);
    
    subplot(1,2,1)
    Y=[T_p05,fliplr(T_p995)];
    Y = Y - s_T_abs(state_ind_T(1));
    hold on
    fill(X,Y,clrmp1(t,:), 'LineWidth', 1);
    scatter(t,s_T_abs(state_ind_T(t)) - s_T_abs(state_ind_T(1)), 'k', 'MarkerFaceColor', 'k') 
    xlabel('time')
    xticks(1:6)
    xticklabels(decades)
    ylabel('degrees C')
    title('Cumulative T Change')
    ylim([-.1 3.75])
    set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' )
    pos1 = pos - [0.15 0 0]; 
    set( yLabelHandle , 'position' , pos1 );
    end
    
    subplot(1,2,2)
    Y=[P_p05,fliplr(P_p995)];
    hold on
    fill(X,Y-P0_abs,clrmp2(t,:), 'LineWidth', 1);
    scatter(t,s_P_abs(state_ind_P(t))-P0_abs, 'k', 'MarkerFaceColor', 'k') 
    xticks(1:6)
    xticklabels(decades)
    xlabel('time')
    ylabel('mm/m')
    title('Cumulative P Change')
    ylim([-18 50])
    set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' )
    pos1 = pos - [0.15 0 0]; 
    set( yLabelHandle , 'position' , pos1 );
    end
        
    frames(t) = getframe(gcf);
    
end

    % Save as movie
    myVideo = VideoWriter('infoOverTime');
    open(myVideo)
    writeVideo(myVideo, frames);
    close(myVideo)

end



end

%% GCM projection plots

cumDeltaT = [zeros(1,21);  cumsum(deltaT(2:end,:))];
cumP = P0_abs *cumprod(deltaP(2:end,:)+1);
cumP = [P0_abs*ones(1,21); cumP];

fig = figure; 
set(fig,'Position', [680 558 1400 750])

t = 1;
    x = t:N+1;
    X=[x,fliplr(x)];
    T_p05 = prctile(T_over_time{t},.01);
    T_p995 = prctile(T_over_time{t},99.9);
    P_p05 = prctile(P_over_time{t},.01);
    P_p995 = prctile(P_over_time{t},99.9);
    
    subplot(1,2,1)
    Y=[T_p05,fliplr(T_p995)];
    Y = Y - s_T_abs(state_ind_T(1));
    hold on
    fill(X,Y,clrmp1(3,:), 'LineWidth', 1);
    scatter(t,s_T_abs(state_ind_T(t)) - s_T_abs(state_ind_T(1)), 'k', 'MarkerFaceColor', 'k') 
    xlabel('time')
    xticks(1:6)
    xticklabels(decades)
    ylabel('degrees C')
    title('Cumulative T Change')
    ylim([-.1 4.5])
    set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' )
    pos1 = pos - [0.15 0 0]; 
    set( yLabelHandle , 'position' , pos1 );
    end
    
    subplot(1,2,2)
    Y=[P_p05,fliplr(P_p995)];
    hold on
    fill(X,Y-P0_abs,clrmp2(3,:), 'LineWidth', 1);
    scatter(t,s_P_abs(state_ind_P(t))-P0_abs, 'k', 'MarkerFaceColor', 'k') 
    xticks(1:6)
    xticklabels(decades)
    xlabel('time')
    ylabel('mm/m')
    title('Cumulative P Change')
    ylim([-18 50])
    set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' )
    pos1 = pos - [0.15 0 0]; 
    set( yLabelHandle , 'position' , pos1 );
    end

subplot(1,2,1)
plot(cumDeltaT, 'k')
hold on
xticks(1:6)
xticklabels(decades)

subplot(1,2,2)
plot(cumP-P0_abs, 'k')
hold on
xticks(1:6)
xticklabels(decades)

fig = figure; 
set(fig,'Position', [680 558 1400 750])
subplot(1,2,1)
plot(cumDeltaT, 'k')
hold on
xlabel('time')
ylim([-.1 4.5])
xticks(1:6)
xticklabels(decades)
ylabel('degrees C')
title('Cumulative T Change')
yLabelHandle = get( gca ,'YLabel' );
pos  = get( yLabelHandle , 'position' )
pos1 = pos - [0.15 0 0]; 
set( yLabelHandle , 'position' , pos1 );
    
subplot(1,2,2)
plot(cumP-P0_abs, 'k')
xlabel('time')
ylabel('mm/m')
title('Cumulative P Change')
ylim([-18 50])
hold on
xticks(1:6)
xticklabels(decades)
yLabelHandle = get( gca ,'YLabel' );
pos  = get( yLabelHandle , 'position' )
pos1 = pos - [0.15 0 0]; 
set( yLabelHandle , 'position' , pos1 );

%% Lambda plot

if false
% Lambda is models X numSamples X decades X states

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

