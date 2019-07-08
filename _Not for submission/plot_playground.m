%Setup
% clear all; close all
load('BMA_results_deltas_2019-01-02.mat')

N = 5;

% Percent change in precip from one time period to next
climParam.P_min = -.3;
climParam.P_max = .3;
climParam.P_delta = .02; 
s_P = climParam.P_min : climParam.P_delta : climParam.P_max;
climParam.P0 = s_P(15);
climParam.P0_abs = 77; %mm/month
M_P = length(s_P);

% Change in temperature from one time period to next
climParam.T_min = 0;
climParam.T_max = 1.5;
climParam.T_delta = 0.05; % deg C
s_T = climParam.T_min: climParam.T_delta: climParam.T_max;
climParam.T0 = s_T(1);
climParam.T0_abs = 26;
M_T = length(s_T);

% Absolute temperature values
T_abs_max = max(s_T) * N;
s_T_abs = climParam.T0_abs : climParam.T_delta : climParam.T0_abs+ T_abs_max;
M_T_abs = length(s_T_abs);
T_bins = [s_T_abs-climParam.T_delta/2 s_T_abs(end)+climParam.T_delta/2];

% Absolute percip values
P_abs_max = max(s_P) * N;
s_P_abs = 66:1:97;
M_P_abs = length(s_P_abs);
P_bins = [s_P_abs-climParam.P_delta/2 s_P_abs(end)+climParam.P_delta/2];

climParam = struct;
climParam.numSamp_delta2abs = 100000;
climParam.numSampTS = 100;
climParam.checkBins = false;

% Percent change in precip from one time period to next
climParam.P_min = -.3;
climParam.P_max = .3;
climParam.P_delta = .02; 
s_P = climParam.P_min : climParam.P_delta : climParam.P_max;
climParam.P0 = s_P(15);
climParam.P0_abs = 77; %mm/month
M_P = length(s_P);

% Change in temperature from one time period to next
climParam.T_min = 0;
climParam.T_max = 1.5;
climParam.T_delta = 0.05; % deg C
s_T = climParam.T_min: climParam.T_delta: climParam.T_max;
climParam.T0 = s_T(1);
climParam.T0_abs = 26;
M_T = length(s_T);

[T_Temp, T_Precip, ~, ~, ~, ~] = bma2TransMat( NUT, NUP, s_T, s_P, N, climParam);

numSamp = 25000;
decades = { '1990', '2010', '2030', '2050', '2070', '2090'};

% Starting point
T0 = s_T(1);
T0_abs = 26;
P0 = s_P(15);
P0_abs = 75;

%% Sample expected transitions over deltas

T_step = s_T(2) - s_T(1);
T_bins = [s_T-T_step/2 s_T(end)+T_step/2];
M_T = length(s_T);

P_step = s_P(2) - s_P(1);
P_bins = [s_P-P_step/2 s_P(end)+P_step/2];
M_P = length(s_P);

T_Temp_del = zeros(M_T, M_T, N);
for t = 1:N
    for index_s_T = 1:M_T
        T_Temp_del(:,index_s_T,t) =  histcounts(NUT(:,t,index_s_T), T_bins, 'Normalization', 'Probability');
    end
end

T_Precip_del = zeros(M_P, M_P, N);
for t = 1:N
    for index_s_P = 1:M_P
        T_Precip_del(:,index_s_P,t) =  histcounts(NUP(:,t,index_s_P), P_bins, 'Normalization', 'Probability');
    end
end

if false
t = 5;
for index_s_T = 1:M_T
    E = sum(T_Temp_del(:,index_s_T,t) .* s_T');
    Tval = s_T(index_s_T);
    Ets(index_s_T) = E;
    Tvalts(index_s_T) = Tval;
    sprintf('E( delT(2050) | delT(2030) = %f ): %f', [Tval E ]) 
end
figure; plot(Tvalts, Ets);

for index_s_P = 1:M_P
    E = sum(T_Precip_del(:,index_s_P,t) .* s_P');
    Pval = s_P(index_s_P);
    Ets(index_s_P) = E;
    Pvalts(index_s_P) = Pval;
    sprintf('E( P(2050) | P(2030) = %f ): %f', [Pval E ]) 
end
figure; plot(Pvalts, Ets);
end

%% Sample bottom up T prob

% Calculate P (27.7 in 2050 | 26.9 in 2030) (delta = .8
% P = sum from i to n of p( delta(3) = .8 | delta(2) = 

% 2: 1990 to 2010, 2: 2010 to 2030, 3: 2030 to 2050, 4: 2050 to 2070, 5: 2070 to 2090

% 
% % Loop over all possible starting values (I extended this to go somehwat
% % below assumed starting mean to take into account varaibilty in mu).

s_T_abs = 25:0.05:s_T_abs(end);
M_T_abs = length(s_T_abs);


p_dist = zeros(M_T_abs,1);  % P dist for all possible values in 2050 given T_2030 = 26.9

for l = 1:M_T_abs
    
    T0 = s_T_abs(l);


    % Sample from Mu to get initial T0 obs, then use this to calculate deltas
    % prior to 1990. We will use to calculate unconditional p(delta1).
    

    indT0 = 2;
    
    R = 1000;
   
    
    T0samp = MUT(:,1,indT0);
    T0samp = T0 + T0samp;
    indsamp = randi(1000,R,1);
    T0samp = T0samp(indsamp);
    T0samp = round2x(T0samp, s_T_abs);
    T0sampDelta = T0samp - T0; 


    T_2030 = 26.9;

    p_dist_by_start = zeros(M_T_abs,1); % P dist for all possible values in 2050 given T_2030 = 26.9 and a starting val;

    % Loop over T values at 2050
    for k = 1:M_T_abs

        T_2050 = s_T_abs(k);

        delta_2030_2050 = T_2050-T_2030;
        delta_1990_2030 = T_2030 - T0;
        ind_delta_2030_2050 = find(abs(s_T - delta_2030_2050) < 0.01);
        p_delta_2030_2050 = zeros(1,M_T);
        p_delta_2010_2030 = zeros(1,M_T);
        p_delta_1990_2010 = zeros(1,M_T);

        p = 0; %prob of this value of T2050

        % loop over all possible detals from 2010 to 2030
        for i = 1:M_T

            % for each possible delta value from 2010 to 2030, calculate what the
            % delta from 1990 to 2010 must have been to get the total delta needed
            % from 1990 to 2030 (.9)
            delta_2010_2030 = s_T(i);
            delta_1990_2010 = delta_1990_2030 - delta_2010_2030;

            % get index in delta state variable for each delta value
            ind_delta_2010_2030 = find(abs(s_T - delta_2010_2030) <0.01);
            ind_delta_1990_2010 = find(abs(s_T - delta_1990_2010) <0.01);

            if isempty(ind_delta_1990_2010) || isempty(ind_delta_2010_2030) || isempty(ind_delta_2030_2050)
                continue; % this combination of deltas is infeasible so don't include
                % (eg if delta_2010_2030 > .9 then delta_1990_2010 must be negative, 
                % which is not possible)
            end

            % get delta transition probabilities 
            p_delta_2030_2050(k) = T_Temp_del(ind_delta_2030_2050, ind_delta_2010_2030, 3);
            p_delta_2010_2030(k) = T_Temp_del(ind_delta_2010_2030, ind_delta_1990_2010, 2);

            % for first time period, need unconditional delta, so average over
            % samples of mu from first time period since state space doesn't tell
            % us about delta prior to 1990
            for j = 1:R
                ind_T0samp_delta = find(abs(T0sampDelta(j) - s_T) < 0.01);
                p_delta_1990_2010_samp(j) = T_Temp_del(ind_delta_1990_2010, ind_T0samp_delta, 1) *1/R;
            end
            p_delta_1990_2010(k) = sum(p_delta_1990_2010_samp);

            % Multiply conditional deltas from each time period and sum over all
            % paths to get total probability 
            p = p + ( p_delta_2030_2050(k) * p_delta_2010_2030(k) * p_delta_1990_2010(k) );

        end
        p_dist_by_start(k) = p;
        
    end
    sum(p_dist_by_start);
    p_dist = p_dist + p_dist_by_start;
    
end
sum(p_dist)
s_T_abs(find(p_dist > 0))
%%
if false
load('Data/Mombasa_TandP.mat'); 
for year = 1:200
    YTij(year,:) = mean(Tij(12*(year-1)+1:12*year,:),1);
    YPij(year,:) = mean(Pij(12*(year-1)+1:12*year,:),1);
end

%% Updating over time with runoff

load('runoff_by_state_Mar16_knnboot_1t')
load('shortageCost_Mar30_damonly_15pen', 'yield', 'unmet_dom', 'shortageCost')


% Example time series 1: wet path

% Set time series
state_ind_P = zeros(1,N);
state_ind_T = zeros(1,N);
state_ind_P(1) =  find(P0_abs==s_P_abs);
state_ind_T(1) = find(T0_abs==s_T_abs);
randGen = true;
state_ind_P(2:N) = [12 17 19 22];
state_ind_T(2:N) = [9 22 35 50];


MAR = cellfun(@(x) mean(mean(x)), runoff);
p = randi(numSamp,N-1);
T_over_time = cell(1,N);
P_over_time = cell(1,N);
MAR_over_time = cell(1,N);

for t = 1:N
    % Sample forward distribution given current state
    T_current = s_T_abs(state_ind_T(t));
    P_current = s_P_abs(state_ind_P(t));
    [T_over_time{t}] = T2forwardSimTemp(T_Temp, s_T_abs, N, t, T_current, numSamp, false);
    [P_over_time{t}] = T2forwardSimTemp(T_Precip, s_P_abs, N, t, P_current, numSamp, false);
    
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
    if randGen
        state_ind_T(t+1) = find(T_over_time{t}(p(t),2)==s_T_abs);
        state_ind_P(t+1) = find(P_over_time{t}(p(t),2)==s_P_abs);
    end
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
    xticks(1:6)
    xticklabels(decades)
    ylabel('degrees C')
    title('Cumulative T Change')
    ylim([-.1 3.75])
    set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' );
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
    ylabel('MCM/y')
    title('Cumulative P Change')
    ylim([-18 20])
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
    ylabel('MCM/y')
    title('Mean Annual Runoff')
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
    ylabel('MCM/y')
    title('Mean Annual Shortage (beyond 10%) ')
    ylim([0 20])
    set(gca,'Units','normalized')
    if t ==1
    yLabelHandle = get( gca ,'YLabel' );
    pos  = get( yLabelHandle , 'position' )
    pos1 = pos - [0.15 0 0]; 
    set( yLabelHandle , 'position' , pos1 );
    end
        
    frames(t) = getframe(gcf);
    
end

end