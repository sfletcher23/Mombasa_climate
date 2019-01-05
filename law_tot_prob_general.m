%Setup
projpath = '/Users/sarahfletcher/Dropbox (MIT)/Mombasa_Climate';
addpath(genpath(projpath))
load('BMA_results_deltas_2019-01-02.mat')

N = 5;

% Delta state variable definition
climParam.T_min = 0;
climParam.T_max = 1.5;
climParam.T_delta = 0.05; % deg C
s_T = climParam.T_min: climParam.T_delta: climParam.T_max;
climParam.T0 = s_T(1);
climParam.T0_abs = 26;
M_T = length(s_T);

% Abs state var def
T_abs_max = max(s_T) * N +climParam.T0_abs;
s_T_abs = 25:0.05:T_abs_max;
M_T_abs = length(s_T_abs);

%% Calculate delta transition probabilities using Nu samples from BMA

T_step = s_T(2) - s_T(1);
T_bins = [s_T-T_step/2 s_T(end)+T_step/2];
M_T = length(s_T);

% P_step = s_P(2) - s_P(1);
% P_bins = [s_P-P_step/2 s_P(end)+P_step/2];
% M_P = length(s_P);

T_Temp_del = zeros(M_T, M_T, N);
for t = 1:N
    for index_s_T = 1:M_T
        T_Temp_del(:,index_s_T,t) =  histcounts(NUT(:,t,index_s_T), T_bins, 'Normalization', 'Probability');
    end
end

% T_Precip_del = zeros(M_P, M_P, N);
% for t = 1:N
%     for index_s_P = 1:M_P
%         T_Precip_del(:,index_s_P,t) =  histcounts(NUP(:,t,index_s_P), P_bins, 'Normalization', 'Probability');
%     end
% end

%% Setup initital time period samples

% Use Mu samples to get initial T0 delta, then use this to calculate deltas
% prior to 1990. We will use to calculate unconditional p(delta1).
    
indT0 = 10;  % This assumes transition prior to 1990; tried several, negligible impact on results
T0samp = MUT(:,1,indT0);
T0samp = round2x(T0samp, s_T);


%% Recursive function to calculate abs transition probabilities 

% function [p_col] = 



%% Calculate sample transition probability using law of total prob

% Calculate P (T=X in 2050 | T=26.9 in 2030)     (Can play with starting val to test)
T_2030 = 26.9; % Assumed in this example


% Time period reference:
% 2: 1990 to 2010, 2: 2010 to 2030, 3: 2030 to 2050, 4: 2050 to 2070, 5: 2070 to 2090



% Loop over all possible starting values (I extended this to go somehwat
% below assumed starting mean to take into account varaibilty in mu).



p_dist = zeros(M_T_abs,1);  % Prob dist for all possible values in 2050 given T_2030 = 26.9

for l = 1:M_T_abs % Starting values in 1990
    
    T0 = s_T_abs(l);
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

        % loop over all possible deltas from 2010 to 2030
        for i = 1:M_T

            % for each possible delta value from 2010 to 2030, calculate what the
            % delta from 1990 to 2010 must have been to get the total delta needed
            % from 1990 to 2030 
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
            R = length(T0samp);
            for j = 1:R
                ind_T0samp= find(abs(T0samp(j) - s_T) < 0.01);
                p_delta_1990_2010_samp(j) = T_Temp_del(ind_delta_1990_2010, ind_T0samp, 1) *1/R;
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
sum(p_dist) % Should ~ equal 1 to be valid prob dist
s_T_abs(find(p_dist > 0)) % Shows 2050 T values with positive probability
