
% Add subfolders to path
if ~isempty(getenv('SLURM_JOB_ID'))
    addpath(genpath('/net/fs02/d2/sfletch/Mombasa_climate'))
else
    addpath(genpath('/Users/sarahfletcher/Documents/MATLAB/Mombasa_Climate'))
end

jobid = getenv('SLURM_JOB_ID');
datetime=datestr(now);
datetime=strrep(datetime,':','_'); %Replace colon with underscore
datetime=strrep(datetime,'-','_');%Replace minus sign with underscore
datetime=strrep(datetime,' ','_');%Replace space with underscore

%% Parameters

runParam = struct;
runParam.N = 5;
runParam.runSDP = true;
runParam.steplen = 20; 
runParam.runRunoff = false;
runParam.runTPts = false;
runParam.runoffPostProcess = false;
runParam.calcTmat = false;
runParam.calcShortage = false;
runParam.runoffLoadName = 'runoff_by_state_comb_Mar2';
runParam.shortageLoadName = 'shortage_costs_28_Feb_2018_17_04_42';
runParam.saveOn = true;

climParam = struct;
climParam.numSamp_delta2abs = 1000;
climParam.numSampTS = 1;
climParam.checkBins = false;

costParam = struct;
costParam.yieldprctl = 80;
costParam.domShortage = 25;
costParam.agShortage = 10;


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
storage = [80 120];

% Actions: Choose dam option in time period 1; expand dam in future time
% periods
a_exp = 0:4; % 0 - do nothing; 1 - build small dam; 2 - build large dam; 3 - build flex dam
            % 4 - expand flex dam
dam_cost = [0 74436346 100192737 74436346*1.08 74436346*.0];

  
%% Calculate climate transition matrix 

% Calculate T_clim

if runParam.calcTmat   
    load('BMA_results_deltap05T_p2P07-Feb-2018 20:18:49.mat')
    [T_Temp, T_Precip, ~, ~, ~, ~] = bma2TransMat( NUT, NUP, s_T, s_P, N, climParam);
    save('T_Temp_Precip', 'T_Temp', 'T_Precip')    
else
    load('T_Temp_Precip') 
end


% Prune state space
for t = 1:N
    index_s_p_time{t} = find(~isnan(T_Precip(1,:,t)));
    index_s_t_time{t} = find(~isnan(T_Temp(1,:,t)));

end


%% Runoff time series for each T,P state

if runParam.runTPts

% Generate T and P timeseries based on T and P mean states
T_ts = cell(M_T_abs,N);
P_ts = cell(M_P_abs,N);
for t = 1:N
    for i = 1:M_T_abs  
        [T_ts{i,t}, ~] = mean2TPtimeseries(t, runParam.steplen, s_P_abs(1), s_T_abs(i), climParam.numSampTS);
    end
    
    for i = 1:M_P_abs  
        [~, P_ts{i,t}] = mean2TPtimeseries(t, runParam.steplen, s_P_abs(i), s_T_abs(1), climParam.numSampTS);
    end
end

end

if runParam.runRunoff 
    
% Generate runoff timeseries - different set for each T,P combination
runoff = cell(M_T_abs, M_P_abs, N);


% Set up parallel 
pc = parcluster('local');
if ~isempty(getenv('SLURM_JOB_ID'))
    parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')));
end

for t = 1:N
    
    % loop over available temp states
    index_s_t_thisPeriod = index_s_t_time{t}; 
    parfor i = 1:length(index_s_t_thisPeriod)
        index_s_t = index_s_t_thisPeriod(i);
        
        runoff_temp = cell(M_P_abs,1);
        
        % loop over available precip states
        index_s_p_thisPeriod = index_s_p_time{t}; 
        for index_s_p = index_s_p_thisPeriod
            
            runoff_temp{index_s_p} = ...
                TP2runoff(T_ts{index_s_t,t}, P_ts{index_s_p,t}, runParam.steplen);

        end
        
        runoff(i, :, t) = runoff_temp;
        
    end
end

savename_runoff = strcat('runoff_by_state_', jobid,'_', datetime);
save(savename_runoff, 'runoff')


if runParam.runoffPostProcess
    % The nature of the parfor loop above saves the runoff timeseries in
    % the wrong temp index; this section corrects that
    runoff_post = cell(M_T_abs, M_P_abs, N);
    for t = 1:N
        
        index_s_p_thisPeriod = index_s_p_time{t}; 
        for index_s_p = index_s_p_thisPeriod
            
            index_s_t_thisPeriod = index_s_t_time{t}; 
            for i= 1:length(index_s_t_thisPeriod)
                
                runoff_post{index_s_t_thisPeriod(i),index_s_p,t} = runoff{i,index_s_p,t};
           
            end
        end
    end
    
    runoff = runoff_post;
    
    savename_runoff = strcat('runoff_by_state_', jobid,'_', datetime);
    save(savename_runoff, 'runoff')
    
end

else
    
load(runParam.runoffLoadName);

end

%% Use reservoir operation model to calculate yield and shortage costs

if runParam.calcShortage

    unmet_ag = zeros(M_T_abs, M_P_abs, length(storage), N);
    unmet_dom = zeros(M_T_abs, M_P_abs, length(storage), N);

    for t = 1:N
        index_s_p_thisPeriod = index_s_p_time{t}; 
        for index_s_p = index_s_p_thisPeriod

            index_s_t_thisPeriod = index_s_t_time{t}; 
            for index_s_t= index_s_t_thisPeriod

                for s = 1:length(storage)

                    [yield_mdl, K, dmd, unmet_dom_mdl, unmet_ag_mdl]  = ...
                        runoff2yield(runoff{index_s_t,index_s_p,t}, T_ts{index_s_t,t}, P_ts{index_s_p,t}, storage(s), runParam, climParam);
                    unmet_ag(index_s_t, index_s_p, s, t) = prctile(sum(unmet_ag_mdl,2),costParam.yieldprctl);
                    unmet_dom(index_s_t, index_s_p, s, t) = prctile(sum(unmet_dom_mdl,2),costParam.yieldprctl);
                end
            end
        end
    end


    shortageCost = (unmet_ag * costParam.agShortage + unmet_dom * costParam.domShortage) * 1E6; 

    savename_shortageCost = strcat('shortage_costs', jobid,'_', datetime);
    save(savename_shortageCost, 'shortageCost')

else
    load(runParam.shortageLoadName);
end
    
%% Backwards Recursion

if runParam.runSDP

% Initialize best value and best action matrices
% Temperature states x precipitaiton states x capacity states, time
V = NaN(M_T_abs, M_P_abs, M_C, N+1);
X = NaN(M_T_abs, M_P_abs, M_C, N+1);

% Terminal period
X(:,:,:,N+1) = zeros(M_T_abs, M_P_abs, M_C, 1);
V(:,:,:,N+1) = zeros(M_T_abs, M_P_abs, M_C, 1);

% Loop over all time periods
for t = linspace(N,1,N)
    
    % Calculate nextV    
    nextV = V(:,:,:,t+1);
          
    % Loop over all states
    
    % Loop over temperature state
    index_s_t_thisPeriod = index_s_t_time{t}; 
    for index_s_t = index_s_t_thisPeriod
        st = s_T_abs(index_s_t);
        
        % Loop over precipitation state
        index_s_p_thisPeriod = index_s_p_time{t}; 
        for index_s_p = index_s_p_thisPeriod
            sp = s_P_abs(index_s_p);
       
            % Loop over capacity expansion state
            for index_s_c = 1:M_C
                sc = s_C(index_s_c);

                bestV = Inf;  % Best value
                bestX = 0;  % Best action

                % Update available actions based on time and whether expansion available
                
                % In first period decide what dam to build
                if t == 1
                    a_exp = 1:3; % build a small, large, or flex dam
                else
                    % In later periods decide whether to expand or not if available
                    switch sc
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

                    % Calculate costs 
                    
                    % Select which capacity is currently available
                    if sc == 1 || sc == 3
                        short_ind = 1;    % small capacity
                    else
                        short_ind = 2; % large capacity
                    end
                    
                    % In first time period, assume have dam built
                    if t == 1
                        if a == 3
                            short_ind = 2; % large capacity
                        else 
                             short_ind = 1;    % small capacity
                        end
                    end
                    
                    cost = shortageCost(index_s_t, index_s_p, short_ind, t) + dam_cost(index_a);
                                      
                   
                    % Calculate transition matrix
                    
                    % Capacity transmat vector
                    T_cap = zeros(1,M_C);
                    if t == 1
                        % In first time period, get whatever dam you picked
                        T_cap(sc) = 1;
                    else
                        % Otherwise, either stay in current or move to expanded
                        switch a
                            case 0
                                T_cap(sc) = 1;
                            case 4
                                T_cap(4) = 1;
                        end                         
                    end

                    % Temperature transmat vector
                    T_Temp_row = T_Temp(:,index_s_t, t)';
                    if sum(isnan(T_Temp_row)) > 0
                        error('Nan in T_Temp_row')
                    end
                    
                    % Precipitation transmat vector
                    T_Precip_row = T_Precip(:,index_s_p, t)';
                    if sum(isnan(T_Precip_row)) > 0
                        error('Nan in T_Precip_row')
                    end
                    
                    % Calculate full transition matrix
                    % T gives probability of next state given current state and actions

                    TRows = cell(3,1);
                    TRows{1} = T_Temp_row;
                    TRows{2} = T_Precip_row;
                    TRows{3} = T_cap;
                    [ T ] = transrow2mat( TRows );

                     % Calculate expected future cost or percentile cost
                    indexNonZeroT = find(T > 0);
                    expV = sum(T(indexNonZeroT) .* nextV(indexNonZeroT));
                    for i = 2:4
                        expV = sum(expV);
                    end

                    stateMsg = strcat('t=', num2str(t), ', st=', num2str(st), ', sp=', num2str(sp), ', sc=', num2str(sc), ', a=', num2str(a))
                    disp(stateMsg)
                    
                   % Check if best decision
                    checkV = cost + expV;
                    if checkV < bestV
                        bestV = checkV;
                        bestX = a;
                    end
                end
            end
            


            % Check that bestV is not Inf
            if bestV == Inf
                error('BestV is Inf, did not pick an action')
            end

            % Save best value and action for current state
            V(index_s_t, index_s_p,index_s_c, t) = bestV;
            X(index_s_t, index_s_p,index_s_c, t) = bestX;

        end
    end
end

if runParam.saveOn
    
    savename_results = strcat('results', jobid,'_', datetime);
    save(savename_results)
    
end


end

    





