
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
runParam.forwardSim = true;
runParam.calcTmat = true;
runParam.calcShortage = true;
runParam.desalOn = true;
runParam.desalCapacity = [60 80];
runParam.runoffLoadName = 'runoff_by_state_Mar16_knnboot_1t';
runParam.shortageLoadName = 'shortage_costs_28_Feb_2018_17_04_42';
runParam.saveOn = true;

climParam = struct;
climParam.numSamp_delta2abs = 100000;
climParam.numSampTS = 20;
climParam.checkBins = false;

costParam = struct;
costParam.yieldprctl = 50;
costParam.domShortage = 15;
costParam.agShortage = 0;
costParam.discountrate = .05;


%% State and Action Definitions 

N = runParam.N;

% Generate state space for climate variables
climParam.P_min = -.3;
climParam.P_max = .3;
climParam.P_delta = .02; %mm/m
s_P = climParam.P_min : climParam.P_delta : climParam.P_max;
climParam.P0 = s_P(15);
climParam.P0_abs = 77;
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
a_exp = 0:4; % 0 - do nothing; 1 - build small; 2 - build large; 3 - build flex
            % 4 - expand flex 
            
infra_cost = zeros(1,length(a_exp));
opex_cost = zeros(1,length(a_exp));
if ~runParam.desalOn
    
    % dam costs
    infra_cost(2) = storage2damcost(storage(1),0);
    infra_cost(3) = storage2damcost(storage(2),0);
    [infra_cost(4), infra_cost(5)] = storage2damcost(storage(1), storage(2));
    percsmalltolarge = (infra_cost(3) - infra_cost(2))/infra_cost(2);
    flexexp = infra_cost(4) + infra_cost(5);
    diffsmalltolarge = infra_cost(3) - infra_cost(2);
    shortagediff = (infra_cost(3) - infra_cost(2))/ (costParam.domShortage * 1e6);
    
else
    % desal capital costs
    [infra_cost(2),~,opex_cost(2)] = capacity2desalcost(runParam.desalCapacity(1),0); % small
    [infra_cost(3),~,opex_cost(3)] = capacity2desalcost(runParam.desalCapacity(2),0); % large
    [infra_cost(4), infra_cost(5), opex_cost(4)] = capacity2desalcost(runParam.desalCapacity(1), runParam.desalCapacity(2));  
    opex_cost(5) = opex_cost(4);

end
infra_cost
opex_cost

fprintf('Large is %.1f MCM shortage \n small is %.1f MCM shoratge \n exp cost is %.1f MCM of shorage \n', ...
    [infra_cost(3) /(costParam.domShortage * 1e6) infra_cost(2) /(costParam.domShortage * 1e6) ...
    infra_cost(5) /(costParam.domShortage * 1e6)])

  
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
    
    
%     % Don't prune
%     index_s_p_time{t} = 1:length(s_P_abs);
%     index_s_t_time{t} = 1:length(s_T_abs);
end


%% Runoff time series for each T,P state

if runParam.runTPts

T_ts = cell(M_T_abs,N);
P_ts = cell(M_P_abs,N);


[Tanom, Panom] = mean2TPtimeseriesMJL_2(1, runParam.steplen, climParam.numSampTS);
for t = 1:N
    
    for i = 1:M_T_abs  
        T_ts{i,t} = Tanom + s_T_abs(i)*ones(size(Tanom));
    end
    
    for i = 1:M_P_abs  
        Ptmp = Panom + s_P_abs(i)*ones(size(Tanom));
        Ptmp(Ptmp<0) = 0;
        P_ts{i,t} = Ptmp;
    end
    
end

savename_runoff = strcat('runoff_by_state_', jobid,'_', datetime);
save(savename_runoff, 'T_ts', 'P_ts')

end

if runParam.runRunoff 
    
% Generate runoff timeseries - different set for each T,P combination
runoff = cell(M_T_abs, M_P_abs, N);


% Set up parallel 
pc = parcluster('local');
if ~isempty(getenv('SLURM_JOB_ID'))
    parpool(pc, str2num(getenv('SLURM_CPUS_ON_NODE')));
end

for t = 1
    
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
save(savename_runoff, 'runoff', 'T_ts', 'P_ts')


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
    save(savename_runoff, 'runoff', 'T_ts', 'P_ts')
    
end

else
    
load(runParam.runoffLoadName);

end

%% Use reservoir operation model to calculate yield and shortage costs

if runParam.calcShortage

    unmet_ag = nan(M_T_abs, M_P_abs, length(storage), N);
    unmet_dom = nan(M_T_abs, M_P_abs, length(storage), N);
    yield = nan(M_T_abs, M_P_abs, length(storage), N);
    desal = cell(M_T_abs, M_P_abs, length(storage));

    for t = 1
        index_s_p_thisPeriod = index_s_p_time{t}; 
        for index_s_p = 1:length(s_P_abs)

            index_s_t_thisPeriod = index_s_t_time{t}; 
            for index_s_t= 1:length(s_T_abs)

                for s = 1:2

                    if ~runParam.desalOn
                        % vary storage
                        [yield_mdl, K, dmd, unmet_dom_mdl, unmet_ag_mdl, desalsupply, desalfill]  = ...
                            runoff2yield(runoff{index_s_t,index_s_p,t}, T_ts{index_s_t,t}, P_ts{index_s_p,t}, storage(s), 0, runParam, climParam);
                    else
                        % large storage, vary desal
                        [yield_mdl, K, dmd, unmet_dom_mdl, unmet_ag_mdl, desalsupply, desalfill]  = ...
                            runoff2yield(runoff{index_s_t,index_s_p,t}, T_ts{index_s_t,t}, P_ts{index_s_p,t}, storage(2), runParam.desalCapacity(s), runParam, climParam);
                    end
                    
                    unmet_dom_90 = max(unmet_dom_mdl - cmpd2mcmpy(186000)*.1, 0);
                    unmet_ag(index_s_t, index_s_p, s, t) = mean(sum(unmet_ag_mdl,2));
                    unmet_dom(index_s_t, index_s_p, s, t) = mean(sum(unmet_dom_90,2));
                    yield(index_s_t, index_s_p, s, t) = mean(sum(yield_mdl,2));
                    desal{index_s_t, index_s_p, s} = desalsupply + desalfill;
                    
                end
            end
        end
    end

    shortageCost =  (unmet_ag * costParam.agShortage + unmet_dom * costParam.domShortage) * 1E6; 
    
    desal_opex = nan(M_T_abs, M_P_abs, length(storage), N);
    for t = 1:N
        discountfactor =  repmat((1+costParam.discountrate) .^ ((t-1)*runParam.steplen+1:t*runParam.steplen), 100, 1);
        desal_opex(:,:,:,t) = cellfun(@(x) mean(opex_cost * x ./ discountfactor, 2), desal);
    end
    

    savename_shortageCost = strcat('shortage_costs', jobid,'_', datetime);
    save(savename_shortageCost, 'shortageCost', 'yield', 'unmet_ag', 'unmet_dom', 'desal_opex')

else
%     load(runParam.shortageLoadName);
end

    
%% Backwards Recursion

if runParam.runSDP

% Initialize best value and best action matrices
% Temperature states x precipitaiton states x capacity states, time
V = NaN(M_T_abs, M_P_abs, M_C, N+1);
X = NaN(M_T_abs, M_P_abs, M_C, N);

% Terminal period
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
                    a_exp_thisPeriod = 1:3; % build a small, large, or flex dam
                else
                    % In later periods decide whether to expand or not if available
                    switch sc
                        case s_C(1) % Small
                            a_exp_thisPeriod = [0];
                        case s_C(2) % Large
                            a_exp_thisPeriod = [0];
                        case s_C(3) % Flex, not expanded
                            a_exp_thisPeriod = [0 4];
                        case s_C(4) % Flex, expanded
                            a_exp_thisPeriod = [0];
                    end
                end
                num_a_exp = length(a_exp_thisPeriod);

                % Loop over expansion action
                for index_a = 1:num_a_exp
                    a = a_exp_thisPeriod(index_a);
                    
                    stateMsg = strcat('t=', num2str(t), ', st=', num2str(st), ', sp=', num2str(sp), ', sc=', num2str(sc), ', a=', num2str(a));
                    disp(stateMsg)

                    % Calculate costs 
                    
                    % Select which capacity is currently available
                    if sc == 1 || sc == 3
                        short_ind = 1;    % small capacity
                    else
                        short_ind = 2; % large capacity
                    end
                    
                    % In first time period, assume have dam built
                    if t == 1
                        if a == 2
                            short_ind = 2; % large capacity
                        else 
                            short_ind = 1;    % small capacity
                        end
                    end
                    
                    % Assume new expansion capacity comes online this period
                    if a == 4
                        short_ind = 2; % large capacity
                    end
                    
                    sCost = shortageCost(index_s_t, index_s_p, short_ind, 1)
                    if t == 1 
                        sCost = 0;  % This is upfront building period
                    end
                  
                    ind_dam = find(a == a_exp);
                    dCost = infra_cost(ind_dam);
                    cost = (sCost + dCost) / (1+costParam.discountrate)^((t-1)*runParam.steplen+1);
                    opex = opex_cost(index_s_t, index_s_p, short_ind, t);
                    cost = cost + opex;
                                      
                   
                    % Calculate transition matrix
                    
                    % Capacity transmat vector
                    T_cap = zeros(1,M_C);
                    if t == 1
                        % In first time period, get whatever dam you picked
                        T_cap(a) = 1;
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
                    
                   % Check if best decision
                    checkV = cost + expV
                    if checkV < bestV
                        bestV = checkV;
                        bestX = a;
                    end
                    
                     % Debug
                    if sc == 3
                        sc
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
end

if runParam.saveOn
    
    savename_results = strcat('results', jobid,'_', datetime);
    save(savename_results)
    
end


end

%% Forward simulation

if runParam.forwardSim
    
% 3 runs: flex, large, small
    
R = 100;
N = runParam.N;

T_state = zeros(R,N);
P_state = zeros(R,N);
C_state = zeros(R,N,4);
action = zeros(R,N,4);
damCostTime = zeros(R,N,4);
shortageCostTime = zeros(R,N,4);
opexCostTime = zeros(R,N,4);
totalCostTime = zeros(R,N,4); 

load('BMA_results_deltap05T_p2P07-Feb-2018 20:18:49.mat', 'MUT', 'MUP')
indT0 = find(s_T_abs == climParam.T0_abs);
indP0 = find(s_P_abs == climParam.P0_abs);
T0samp = MUT(:,1,indT0);
P0samp = MUP(:,1,indP0);
T0samp = climParam.T0_abs + T0samp;
P0samp = exp(P0samp)* climParam.P0_abs;
indsamp = randi(1000,R,1);
T0samp = T0samp(indsamp);
P0samp = P0samp(indsamp);
T0samp = round2x(T0samp, s_T_abs);
P0samp = round2x(P0samp, s_P_abs);

T_state(:,1) = T0samp;
P_state(:,1) = P0samp;
C_state(:,1,1) = 3; % Always flex
C_state(:,1,2) = 2; % Always large
C_state(:,1,3) = 1; % Always small
C_state(:,1,4) = 1; % Choose based on policy

for k = 1:4
    for i = 1:R
        for t = 1:N

            % Choose best action given current state
            index_t = find(T_state(i,t) == s_T_abs);
            index_p = find(P_state(i,t) == s_P_abs);
            index_c = find(C_state(i,t,k) == s_C);
            
            % In flex case follow exp policy, otherwise restrict to large or
            % small and then no exp
            if t==1
                switch k
                    case 1
                        action(i,t,k) = 3;
                    case 2
                        action(i,t,k) = 2;
                    case 3
                        action(i,t,k) = 1;
                    case 4
                        action(i,t,k) =  X(index_t, index_p, index_c, t);
                end
            else 
                switch k
                    case 1
                        action(i,t,k) = X(index_t, index_p, index_c, t);
                    case 2
                        action(i,t,k) = 0;
                    case 3
                        action(i,t,k) = 0;
                    case 4
                        action(i,t,k) =  X(index_t, index_p, index_c, t);
                end
            end
            
            % Save costs of that action

            % Get current capacity and action
            sc = C_state(i,t,k);
            a = action(i,t,k);

            % Select which capacity is currently available
            if sc == 1 || sc == 3
                short_ind = 1;    % small capacity
            else
                short_ind = 2; % large capacity
            end

            % In first time period, assume have dam built
            if t == 1
                if a == 2
                    short_ind = 2; % large capacity
                else 
                    short_ind = 1;    % small capacity
                end
            end

            % Assume new expansion capacity comes online this period
            if a == 4
                short_ind = 2; % large capacity
            end

            % Get shortage and dam costs
            shortageCostTime(i,t,k) = shortageCost(index_t, index_p, short_ind, 1);
            if t == 1 
                shortageCostTime(i,t,k) = 0;  % This is upfront building period
            end
            ind_dam = find(a == a_exp);
            damCostTime(i,t,k) = infra_cost(ind_dam);
            opexCostTime(i,t,k) = opex_cost(index_t, index_p, short_ind, t);
            totalCostTime(i,t,k) = (shortageCostTime(i,t,k) + damCostTime(i,t,k)) / (1+costParam.discountrate)^((t-1)*runParam.steplen+1);
            totalCostTime(i,t,k) = totalCostTime(i,t,k) + opexCostTime(i,t,k);
            

                    cost = (sCost + dCost) / (1+costParam.discountrate)^((t-1)*runParam.steplen+1);
                    opex = opex_cost(index_s_t, index_s_p, short_ind, 1);
                    cost = cost + opex;

            % Simulate transition to next state
            % Capacity transmat vector
            T_cap = zeros(1,M_C);
            if t == 1
                % In first time period, get whatever dam you picked
                T_cap(a) = 1;
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
            T_Temp_row = T_Temp(:,index_t, t)';
            if sum(isnan(T_Temp_row)) > 0
                error('Nan in T_Temp_row')
            end

            % Precipitation transmat vector
            T_Precip_row = T_Precip(:,index_p, t)';
            if sum(isnan(T_Precip_row)) > 0
                error('Nan in T_Precip_row')
            end
            % Combine trans vectors into matrix
            TRows = cell(3,1);
            TRows{1} = T_Temp_row;
            TRows{2} = T_Precip_row;
            TRows{3} = T_cap;
            [ T_current ] = transrow2mat( TRows );

            % Simulate next state
            if t < N
                T_current_1D = reshape(T_current,[1 numel(T_current)]);
                T_current_1D_cumsum = cumsum(T_current_1D);
                p = rand();
                index = find(p < T_current_1D_cumsum,1);
                [ind_s1, ind_s2, ind_s3] = ind2sub(size(T_current),index);
                    % Test sample
                    margin = 1e-10;
                    if (T_current(ind_s1, ind_s2, ind_s3) < margin)
                        error('Invalid sample from T_current')
                    end
                
                if k == 1
                    T_state(i,t+1,k) = s_T_abs(ind_s1);
                    P_state(i,t+1,k) = s_P_abs(ind_s2);
                end
                C_state(i,t+1,k) = s_C(ind_s3);

            end



        end

    end
end


% if runParam.saveOn
%     
%     savename_results = strcat('results', jobid,'_', datetime);
%     save(savename_results)
%     
% end



end





