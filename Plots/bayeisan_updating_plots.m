%Setup
load('BMA_results_deltap05T_p2P07-Feb-2018 20:18:49.mat')


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

load('IPCC_CRU')

%% Sample through T_Temp to get delta time series
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
P_delta_over_time = s_P(state_ind_P);

%% Calculate absolute values from deltas

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

%% Updating over time

p = randi(numSamp,N-1);
state_ind_P = zeros(1,N);
state_ind_T = zeros(1,N);
state_ind_P(1) =  find(P0_abs==s_P_abs);
state_ind_T(1) = find(T0_abs==s_T_abs);
% T0_ind = randi(31,1,1);
% state_ind_T(1) = T0_ind;
T_over_time = cell(1,N);
P_over_time = cell(1,N);
MAR_over_time = cell(1,N);

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

fig = figure;
[clrmp1]=cbrewer('seq', 'Reds', N);
[clrmp2]=cbrewer('seq', 'Blues', N);
set(fig,'Position', [680 558 1400 750])


for t =1:N
    x = t:N+1;
    X=[x,fliplr(x)];
    T_p01 = prctile(T_over_time{t},.01);
    T_p995 = prctile(T_over_time{t},99.9);
    P_p01 = prctile(P_over_time{t},.01);
    P_p995 = prctile(P_over_time{t},99.9);
    
    subplot(1,2,1)
    Y=[T_p01,fliplr(T_p995)];
    Y = Y - s_T_abs(state_ind_T(1));
    hold on
    fill(X,Y,clrmp1(t,:), 'LineWidth', 1);
    scatter(t,s_T_abs(state_ind_T(t)) - s_T_abs(state_ind_T(1)), 'k', 'MarkerFaceColor', 'k') 
    xlabel('time')
    set(gca,'XTick', [1:6], 'XTickLabels',decades)
    %xticks(1:6)
    %xticklabels(decades)
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
    Y=[P_p01,fliplr(P_p995)];
    hold on
    fill(X,Y-P0_abs,clrmp2(t,:), 'LineWidth', 1);
    scatter(t,s_P_abs(state_ind_P(t))-P0_abs, 'k', 'MarkerFaceColor', 'k') 
    set(gca,'XTick', [1:6], 'XTickLabels',decades)
    %xticks(1:6)
    %xticklabels(decades)
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

if false
    % Save as movie
    myVideo = VideoWriter('infoOverTime');
    open(myVideo)
    writeVideo(myVideo, frames);
    close(myVideo)

end

%% Updating over time with runoff

load('runoff_by_state_Mar16_knnboot_1t')
load('shortageCost_Mar30_damonly_15pen', 'yield', 'unmet_dom', 'shortageCost')
MAR = cellfun(@(x) mean(mean(x)), runoff);



p = randi(numSamp,N-1);
state_ind_P = zeros(1,N);
state_ind_T = zeros(1,N);
state_ind_P(1) =  find(P0_abs==s_P_abs);
state_ind_T(1) = find(T0_abs==s_T_abs);
% T0_ind = randi(31,1,1);
% state_ind_T(1) = T0_ind;
T_over_time = cell(1,N);
P_over_time = cell(1,N);
MAR_over_time = cell(1,N);

for t = 1:N
    % Sample forward distribution given current state
    T_current = s_T_abs(state_ind_T(t));
    P_current = s_P_abs(state_ind_P(t));
    [T_over_time{t}] = T2forwardSimTemp(T_Temp_abs, s_T_abs, N, t, T_current, numSamp, false);
    [P_over_time{t}] = T2forwardSimTemp(T_Precip_abs, s_P_abs, N, t, P_current, numSamp, false);
    
    % Lookup MAR and yield for forward distribution
    T_ind = arrayfun(@(x) find(x == s_T_abs), T_over_time{t});
    P_ind = arrayfun(@(x) find(x == s_P_abs), P_over_time{t});
    [~,t_steps] = size(T_ind);
    MAR_over_time{t} = zeros(size(T_ind));
    yield_over_time{t} = zeros(size(T_ind));
    for i = 1:numSamp
        for j = 1:t_steps   
            MAR_over_time{t}(i,j) = MAR(T_ind(i,j), P_ind(i,j), 1);
            yield_over_time{t}(i,j) = unmet_dom(T_ind(i,j), P_ind(i,j),1, 1) ;   % 80 MCM storage
        end
    end
    
    % Sample next time period
    state_ind_T(t+1) = find(T_over_time{t}(p(t),2)==s_T_abs);
    state_ind_P(t+1) = find(P_over_time{t}(p(t),2)==s_P_abs);
end

fig = figure;
[clrmp1]=cbrewer('seq', 'Reds', N);
[clrmp2]=cbrewer('seq', 'Blues', N);
[clrmp3]=cbrewer('seq', 'Greens', N);
[clrmp4]=cbrewer('seq', 'Purples', N);
set(fig,'Position', [680 558 1400 750])


for t =1:N
    x = t:N+1;
    X=[x,fliplr(x)];
    T_p01 = prctile(T_over_time{t},.01);
    T_p995 = prctile(T_over_time{t},99.9);
    P_p01 = prctile(P_over_time{t},.01);
    P_p995 = prctile(P_over_time{t},99.9);
    MAR_p01 = prctile(MAR_over_time{t},.01);
    MAR_p995 = prctile(MAR_over_time{t},99.9);
    yield_p01 = prctile(yield_over_time{t},.01);
    yield_p995 = prctile(yield_over_time{t},99.9);
    
    subplot(2,2,1)
    Y=[T_p01,fliplr(T_p995)];
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
    
    subplot(2,2,2)
    Y=[P_p01,fliplr(P_p995)];
    hold on
    fill(X,Y-P0_abs,clrmp2(t,:), 'LineWidth', 1);
    scatter(t,s_P_abs(state_ind_P(t))-P0_abs, 'k', 'MarkerFaceColor', 'k') 
    xticks(1:6)
    xticklabels(decades)
    xlabel('time')
    ylabel('MCM/y')
    title('Cumulative P Change')
    ylim([-18 50])
    set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' )
    pos1 = pos - [0.15 0 0]; 
    set( yLabelHandle , 'position' , pos1 );
    end
    
    subplot(2,2,3)
    Y=[MAR_p01,fliplr(MAR_p995)];
    hold on
    fill(X,Y,clrmp3(t,:), 'LineWidth', 1);
    scatter(t,MAR_over_time{t}(1,1), 'k', 'MarkerFaceColor', 'k') 
    xticks(1:6)
    xticklabels(decades)
    xlabel('time')
    ylabel('MCM/y')
    title('Mean Annual Runoff')
%     ylim([-18 50])
    set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' )
    pos1 = pos - [0.15 0 0]; 
    set( yLabelHandle , 'position' , pos1 );
    end
    
    subplot(2,2,4)
    Y=[yield_p01,fliplr(yield_p995)];
    hold on
    fill(X,Y,clrmp4(t,:), 'LineWidth', 1);
    scatter(t,yield_over_time{t}(1,1), 'k', 'MarkerFaceColor', 'k') 
    xticks(1:6)
    xticklabels(decades)
    xlabel('time')
    ylabel('MCM/y')
    title('Mean Annual Shortage (beyond 10%) ')
%     ylim([-18 50])
    set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' )
    pos1 = pos - [0.15 0 0]; 
    set( yLabelHandle , 'position' , pos1 );
    end
        
    frames(t) = getframe(gcf);
    
end



set(findall(fig.Children,'Type', 'Scatter'), 'SizeData', 10)
for i = 1:4
    ax = fig.Children(i);
    ax.YLabel.Position = ax.YLabel.Position - [.18 0 0];
end

for i = 1:5
    ax = fig.Children(i);
    ax.XLabel.Position = ax.XLabel.Position - [0 .18 0];
end

%% Add individual GCM projections as deltas

load('GCMdata_deltaTandP')

cumDeltaT = [zeros(1,21);  cumsum(deltaT(2:end,:))];
cumP = P0_abs *cumprod(deltaP(2:end,:)+1);
cumP = [P0_abs*ones(1,21); cumP];



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


fig1 = figure; 
set(fig1,'Position', [680 558 1400 750])

t = 1;
    x = t:N+1;
    X=[x,fliplr(x)];
    T_p01 = prctile(T_over_time{t},.01);
    T_p995 = prctile(T_over_time{t},99.9);
    P_p01 = prctile(P_over_time{t},.01);
    P_p995 = prctile(P_over_time{t},99.9);
    
    subplot(2,2,1)
    Y=[T_p01,fliplr(T_p995)];
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
    
    subplot(2,2,2)
    Y=[P_p01,fliplr(P_p995)];
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

subplot(2,2,1)
plot(cumDeltaT, 'k')
hold on
xticks(1:6)
xticklabels(decades)

subplot(2,2,2)
plot(cumP-P0_abs, 'k')
hold on
xticks(1:6)
xticklabels(decades)

%%
% fig1 = figure; 
subplot(2,2,3) 
boundedline([2005:1:2100]',IPCC.meanT(105:200)',[1.64*IPCC.stdT(105:200)' 1.64*IPCC.stdT(105:200)'],'alpha','cmap',[0.7,0,0]); 
hold on; 
boundedline([1911:1:2005]',IPCC.meanT(11:105)',[1.64*IPCC.stdT(11:105)' 1.64*IPCC.stdT(11:105)'],'alpha','cmap',[0.1,0.1,0.1]); 
hold on; 
plot([1911:1:2005]',CRU.Tsmooth(11:105)-mean(CRU.T(86:105)),'k','LineWidth',2)
xlim([1911,2100]); box on; h = vline(2005,'k'); 
set(gca,'FontSize', 6); xlabel('Year'); ylabel('\Delta T (relative to 1986-2005)');

subplot(2,2,4) 
boundedline([2005:1:2100]',IPCC.meanP(105:200)',[1.64*IPCC.stdP(105:200)' 1.64*IPCC.stdP(105:200)'],'alpha','cmap',[0,0,0.7]); 
hold on; 
boundedline([1911:1:2005]',IPCC.meanP(11:105)',[1.64*IPCC.stdP(11:105)' 1.64*IPCC.stdP(11:105)'],'alpha','cmap',[0.1,0.1,0.1]); 
hold on; 
plot([1911:1:2005]',CRU.Psmooth(11:105)-mean(CRU.P(86:105)),'k','LineWidth',2)
xlim([1911,2100]); box on; h = vline(2005,'k'); ylabel('\Delta P (relative to 1986-2005)');xlabel('Year')
set(gca,'FontSize', 6)
figure_width = 7; % in inches
figure_height = 3; % in inches
font_size = 6; % in points (relative to figure size in inches)
savingind = 1;
figname = 'IPCCbounded';
% figurefun(figure_width, figure_height,font_size,savingind,figname,fig1)