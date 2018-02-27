
% Add subfolders to path
if ~isempty(getenv('SLURM_JOB_ID'))
    addpath(genpath('/net/fs02/d2/sfletch/Mombasa_climate'))
else
    addpath(genpath('/Users/sarahfletcher/Documents/MATLAB/Mombasa_Climate'))
end


%% Parameters

runParam = struct;
runParam.N = 5;
runParam.runSDP = false;
runParam.steplen = 20; 
runParam.runRunoff = true;

climParam = struct;
climParam.numSamp_delta2abs = 10000;
climParam.numSampTS = 1000;
climParam.checkBins = false;


%% State and Action Definitions 

N = runParam.N;

% Generate state space for climate variables
climParam.P_min = -.3;
climParam.P_max = .3;
climParam.P_delta = .02; %mm/m
s_P = climParam.P_min : climParam.P_delta : climParam.P_max;
climParam.P0 = s_P(15);
climParam.P0_abs = 75;
M_P = length(s_P);

climParam.T_min = 0;
climParam.T_max = 1.5;
climParam.T_delta = 0.05; % deg C
s_T = climParam.T_min: climParam.T_delta: climParam.T_max;
climParam.T0 = s_T(1);
climParam.T0_abs = 26;
M_T = length(s_T);


T_abs_max = max(s_T) * N;
s_T_abs = climParam.T0_abs : climParam.T_delta : climParam.T0_abs+ T_abs_max;
M_T_abs = length(s_T_abs);
T_bins = [s_T_abs-climParam.T_delta/2 s_T_abs(end)+climParam.T_delta/2];
T_Temp_abs = zeros(M_T_abs,M_T_abs,N);


P_abs_max = max(s_P) * N;
s_P_abs = 66:1:97;
M_P_abs = length(s_P_abs);
P_bins = [s_P_abs-climParam.P_delta/2 s_P_abs(end)+climParam.P_delta/2];
T_Precip_abs = zeros(M_P_abs,M_P_abs,N);


% State space for capacity variables
s_C = 1:4; % 1 - small;  2 - large; 3 - flex, no exp; 4 - flex, exp
M_C = length(s_C);

% Actions: Choose dam option in time period 1; expand dam in future time
% periods
a_exp = 0:4; % 0 - do nothing; 1 - build small dam; 2 - build large dam; 3 - build flex dam
            % 4 - expand flex dam

  
%% Calculate climate transition matrix 


% Calculate T_clim

load('BMA_results_deltap05T_p2P07-Feb-2018 20:18:49.mat')

[T_Temp, T_Precip, ~, ~, ~, ~] = bma2TransMat( NUT, NUP, s_T, s_P, N, climParam);


% Prune state space
for t = 1:N
    index_s_p_time{t} = find(~isnan(T_Precip(1,:,t)));
    index_s_t_time{t} = find(~isnan(T_Temp(1,:,t)));

end


%% Runoff time series for each T,P state

if runParam.runRunoff 

% Generate T and P timeseries based on T and P mean states
T_ts = cell(M_T,N);
P_ts = cell(M_P,N);
for t = 1:N
    for i = 1:M_T  % conveniently, there are the same number of states for T and P right now
        [T_ts{i,t}, P_ts{i,t}] = mean2TPtimeseries(t, runParam.steplen, s_P_abs(i), s_T_abs(i), climParam.numSampTS);
    end
end

% Generate runoff timeseries - different set for each T,P combination
runoff = cell(M_T, M_P, N);

for t = 1:N
    
    % loop over available temp states
    index_s_t_thisPeriod = index_s_t_time{t}; 
    for index_s_t = index_s_t_thisPeriod
        st = s_T(index_s_t);
        
        % loop over available precip states
        index_s_p_thisPeriod = index_s_p_time{t}; 
        for index_s_p = index_s_p_thisPeriod
            sp = s_P(index_s_p); 
            
            runoff{index_s_t, index_s_p, t} = ...
                TP2runoff(T_ts{index_s_t,t}, P_ts{index_s_p,t}, runParam.steplen);

        end
    end
end

save('runoff_by_state', 'runoff')

end

%% Backwards Recursion

if runParam.runSDP

% Initialize best value and best action matrices
% Temperature states x precipitaiton states x capacity states, time
V = NaN(M_T, M_P, M_C, N+1);
X = NaN(M_T, M_P, M_C, N+1);

% Terminal period
X(:,:,:,N+1) = zeros(M_T, M_P, M_C, 1);
V(:,:,:,N+1) = zeros(M_T, M_P, M_C, 1);

% Loop over all time periods
for t = linspace(N,1,N)
    
    % Calculate nextV    
    nextV = V(:,:,t+1);
          
    % Loop over all states
    
    % Loop over temperature state
    index_s_t_thisPeriod = index_s_t_time{t}; 
    for index_s_t = index_s_t_thisPeriod
        st = s_T(index_s_t);
        
        % Loop over precipitation state
        index_s_p_thisPeriod = index_s_p_time{t}; 
        for index_s_p = index_s_p_thisPeriod
            sp = s_P(index_s_p);
       
            % Loop over capacity expansion state
            for index_s_c = 1:M_C
                sc = s_expand(index_s2);

                bestV = Inf;  % Best value
                bestX = 0;  % Best action

                % Update available actions based on time and whether expansion available
                
                % In first period decide what dam to build
                if t == 1
                    a_exp = 1:3; % build a small, large, or flex dam
                else
                    % In later periods decide whether to expand or not if available
                    switch s2
                        case s_C(1) % Small
                            a_exp = [0];
                        case s_C(2) % Large
                            a_exp = [0];
                        case s_C(3) % Flex, not expanded
                            a_exp = [0 4];
                        case s_C(4) % Flex, expanded
                            a_exp = [0];
                    end
                end
                num_a_exp = length(a_exp);

                % Loop over expansion action
                for index_a = 1:num_a_exp
                    a = a_exp(index_a);

                    % Generate 
                    
                    
                    
                    % Updated to here so far!!!
                    
                    
                    
                    
                    
                    % Calculate demand
%                     demandThisPeriod = demand(water, population(t), t, gwParam);
                    demandThisPeriod = gwParam.pumpingRate;

                    % Calculate cost and shortages this period
                    [ cost, ~,~, ~,~, ~, ~, ~, ~] = cost_supply_func( a1, a, s1, s2, costParam, water, gwParam, t, demandThisPeriod, runParam.capacityDelay, exp_vectors, false);

                    % Calculate transition matrix
                    
                    % If stop pumping, move to state -1. Otherwise, use
                    % T_gw calculated above. 
                    
                    switch a
                        case 0
                            T_gw = zeros(1,gw_M);
                            T_gw(1) = 1;
                        case 1
                            T_gw = T_gw_all(:,index_s1,t)';
                    end

                    % Get transmat vector for next expansion state
                    % (deterministic)                  
                    T_expand = zeros(1,exp_M);
                    
                    if runParam.capacityDelay
                        T_exp_online_ind = subindex_s2(1);
                        T_exp_delay1_ind = subindex_s2(2);
                        T_exp_delay2_ind = subindex_s2(3);
                        % Move delayed capacity to online
                        if subindex_s2(3) == 2  % Move big plant from delay2 to delay 1
                            T_exp_delay1_ind = 4;
                            T_exp_delay2_ind = 1;
                        elseif subindex_s2(2) > 1
                            T_exp_online_ind = T_exp_online_ind + (subindex_s2(2) - 1);
                            T_exp_delay1_ind = 1;
                        end
                        % Add new capacity to delay
                        if a == 1
                            T_exp_delay1_ind = 2;
                        elseif a == 2
                            T_exp_delay2_ind = 2;
                        end
                        temp_index = vectorIndex([T_exp_online_ind T_exp_delay1_ind T_exp_delay2_ind], {s_exp_on', s_exp_delay1', s_exp_delay2'});
                        exp_index = find(s_expand == temp_index);
                        T_expand(exp_index) = 1;
                    else
                        if a == 0
                            T_expand(index_s2) = 1; % Stay in current state
                        elseif a == 1
                            T_expand(index_s2 + 1) = 1; % Move up one state
                        elseif a == 2
                            T_expand(index_s2 + 3) = 1; % Move up three states
                        end
                    end
                    
                    % Calculate full transition matrix
                    % T gives probability of next state given
                    % current state and actions

                    TRows = cell(2,1);
                    TRows{1} = T_gw;
                    TRows{2} = T_expand;
                    [ T ] = transrow2mat( TRows );

                     % Calculate expected future cost or percentile cost
                    indexNonZeroT = find(T > 0);
                    if runParam.percentile
                        nonZeroV = nextV(indexNonZeroT);
                        [sortedNextV, indexSort] = sort(reshape(nonZeroV, 1, []));
                        TnonZero = T(indexNonZeroT);
                        indexPrcnV = find(cumsum(TnonZero(indexSort)) > runParam.percentile/100, 1);
                        expV = nonZeroV(indexPrcnV);
                    else
                        expV = sum(T(indexNonZeroT) .* nextV(indexNonZeroT));
                        for i = 2:4
                            expV = sum(expV);
                        end
                    end
                    
                    stateMsg = strcat('t=', num2str(t), ', s1=', num2str(s1), ', a1=', num2str(a1), ', s2=', num2str(s2), ', a2=', num2str(a))
                    disp(stateMsg)
                   % Check if best decision
                    checkV = cost + expV;
                    if checkV < bestV
                        bestV = checkV;
                        bestX1 = a1;
                        bestX2 = a;
                    end
                end
            end
            


            % Check that bestV is not Inf
            if bestV == Inf
                error('BestV is Inf, did not pick an action')
            end

            % Save best value and action for current state
            V(index_s1, index_s2, t) = bestV;
            X1(index_s1, index_s2, t) = bestX1;
            X2(index_s1, index_s2, t) = bestX2;

        end
    end
end


end

    





