function [ V, X1, X2, T_gw_all, cumTgw, s_gw, s_expand, exp_vectors, lowestCost, lowestCostAction ] = ...
    sdp_gw( runParam, costParam, gwParam, water, datetime )

% Run SDP for groundwater model. 
V = [];
X1 = [];
X2 = [];

%% State and Action Definitions for Groundwater 

% Generate state space for groundwater head and demand range
[s_gw, gw_M] = gen_water_growth_states(gwParam);

% Actions: Stop pumping groundwater (0), continue pumping (1)
a_gw_available = [0 1];
a_gw_unavailable = [0];
  

%% Desalination Expansions: State Definitions and Actions

% a2 desal actions: 0 no expand, 1 expand small, 2 expand large
a_expand_available = [0 1 2];
a_expand_unavailable = 0;

% State definition: volume of additional capacity
maxNumSmallExp = 3;
maxNumLargeExp = 1;

% Check that large capacity is a multiple of small capacity
if mod(water.desal_capacity_expansion.large , water.desal_capacity_expansion.small) ~= 0
    error('Large capacity is not a multiple of small capacity')
end

% Get max capacity, state space between 0 and max cap in steps of small capacity
maxExpCap = water.desal_capacity_expansion.large;
s_expand = 0:water.desal_capacity_expansion.small:maxExpCap;
exp_M = length(s_expand); % Desalination expanded = 2

% Add capacity delay to state space
exp_vectors = [];
if runParam.capacityDelay
   s_exp_on = s_expand;
   s_exp_delay1 = s_expand;
   s_exp_delay2 = [s_expand(1) s_expand(2)];
   % Index for feasible expansion state combinations - max total 3v across substates
    s_expand = [1:7 9 10 13 17];
    exp_M = length(s_expand);
    exp_vectors = {s_exp_on' s_exp_delay1' s_exp_delay2'};
end

% Set cost function
if runParam.oldCost
    cost_supply_func = str2func('supplyAndCost_old');
else
    cost_supply_func = str2func('supplyAndCost');
end

%% Get K and S samples and use to prune state space

N = runParam.N;



    
% Load parameter and drawdown samples
load('T_gw_inputs_Dec4','drawdown')

% Get min and max hydrograph
drawdown_max = [ 71.8024349691912,	98.9996385259814,	121.941450801024,	141.343687632105,	157.856744357095,	172.048031095633,	184.398168603787, ...
                     195.305381025659,	205.094059319596,	214.024853008525,	222.304709021687,	230.096011944031,	237.524453528662,	244.685548197193, ...
                     251.649884766433,	258.467318232683,	265.170400521298,	271.777451521216,	278.295783299501,	284.725670246833,	291.065599430146, ...
                     297.318968271482,	303.501574870096,	309.648071618181,	315.814627673233,	322.075373217118,	328.512369315768,	335.202000802936, ...
                     342.202671675442,	349.547910713690];
drawdown_min = [7.90152696173453,	10.8863953854826,	13.3341292794019,	15.3663710702284,	17.0783085592299,	18.5440767178084,	19.8211011075899, ...
                    20.9535357329021,	21.9749402601961,	22.9103236238453,	23.7776598095863,	24.5889611257012,	25.3509773897663,	26.0655789803655, ...
                    26.7298813924321,	27.3361843409086,	27.8718371562427,	28.3192122981390,	28.6560732294635,	28.8567461499751,	28.8945904598435, ...
                    28.7461885505084,	28.3972642466096,	27.8494647933049,	27.1259694602479,	26.2731049674530,	25.3557618531984,	24.4468338880858, ...
                    23.6140170959438,	22.9088843642109] ;


for t = 1:N
    indexValidState = ceil(drawdown_min(t)) <= s_gw & s_gw <= floor(drawdown_max(t));
    index_s_gw_time{t} = find(indexValidState);
    index_s_gw_time{t} =  [1 index_s_gw_time{t}];
end

if runParam.calculateTgw
% Calculate Groudnwater transition matrix when pumping

    % Get transmat vector for gw when pumping for current gw state

    T_gw_all = zeros(gw_M, gw_M, N);

    for t =1:N
        for index_s1 = 1:gw_M
            s1_now = s_gw(index_s1);
            if ismember(index_s1, index_s_gw_time{t})
                dd_input = drawdown{index_s1,t};
                [T_gw_temp] = gw_transrow_numint(gwParam, s1_now, s_gw, dd_input )';
            else
                T_gw_temp = NaN;
            end
            T_gw_all(:,index_s1,t) = T_gw_temp';
        end
    end    
    
    %save(strcat('T_gw_',datetime), 'T_gw_all')
    
else
    load(gwParam.TgwLoadName)
    T_gw_all = T_gw;
end

% Calculate expected total drawdown for each state
cumTgw = zeros(gw_M, N);
T_gw_temp = T_gw_all;
indexNan = isnan(T_gw_temp);
T_gw_temp(indexNan) = 0;
for t = linspace(N,1,N)
    for index_s1 = 1:gw_M
        s1 = s_gw(index_s1);
        if t == N 
            cumTgw(index_s1,t) = T_gw_temp(1,index_s1,t);
        else
            cumTgw(index_s1,t) = sum(T_gw_temp(:,index_s1,t) .* cumTgw(:,t+1));
        end
    end
end


%% Backwards Recursion

if runParam.runSDP

% Initialize best value and best action matrices
% Groundwater states x desal states x time
V = NaN(gw_M, exp_M, N+1);
X1 = NaN(gw_M, exp_M, N+1);
X2 = NaN(gw_M, exp_M, N+1);

% Terminal period
X1(:,:,N+1) = zeros(gw_M, exp_M, 1);
X2(:,:,N+1) = zeros(gw_M, exp_M, 1);
V(:,:,N+1) = zeros(gw_M, exp_M, 1);

% Loop over all time periods
for t = linspace(N,1,N)
    % Calculate nextV    
    nextV = V(:,:,t+1);
          
    % Loop over all states
    
    % Loop over groundwater state: 1 is depleted, M1 is full
    index_s_gw_thisPeriod = index_s_gw_time{t}; 
    for index_s1 = index_s_gw_thisPeriod
        s1 = s_gw(index_s1);
       
        % Loop over expansion state: 1 is unexpanded, 2 is expanded
        for index_s2 = 1:exp_M
            s2 = s_expand(index_s2);
            
            if runParam.capacityDelay
                subindex_s2 = linIndex2VecIndex(s2, {s_exp_on', s_exp_delay1', s_exp_delay2'});
            end

            bestV = Inf;
            bestX1= 0;  % Groundwater action and expansion action
            bestX2= 0;
            
            % Update available actions based on whether gw available
            if s1 == 200
                a_gw = a_gw_unavailable;    % unavailble bc depleted
            elseif s1 == -1
                a_gw = a_gw_unavailable;    % unavailble bc turned off
            else
                a_gw = a_gw_available;
            end

            % Update available actions based on whether expansion available
            if runParam.capacityDelay
                if subindex_s2(1) == 4 || subindex_s2(2) == 4 || subindex_s2(3) == 2 ... % Max capacity in any of subcategories
                        || (subindex_s2(1) + subindex_s2(2)) >= 5  % Max capacity in online + delay 1
                    a_expand = [0]; % If max capacity online or waiting to come online, can't expand
                elseif (subindex_s2(1) + subindex_s2(2)) == 4 ...   % 2 small units online or in delay 1
                        || (subindex_s2(1) + subindex_s2(2)) == 3    % 1 small units online or in delay 1
                    a_expand = [0 1];
                elseif (subindex_s2(1) + subindex_s2(2) + subindex_s2(3)) == 3
                    a_expand = [0 1 2];
                else
                    error('Some combination was not included!')
                end
                
            else
                switch s2
                    case s_expand(1)
                        a_expand = [0 1 2];
                    case s_expand(2)
                        a_expand = [0 1];
                    case s_expand(3)
                        a_expand = [0 1];
                    case s_expand(4)
                        a_expand = [0];
                end
            end

            num_a_gw = length(a_gw);
            num_a_expand = length(a_expand);

            % Loop over all actions
            % Loop over groundwater pumping action
            for index_a1 = 1:num_a_gw
                a1 = a_gw(index_a1);

                % Loop over expansion action: 1 is do not expand, 2 is expand
                for index_a2 = 1:num_a_expand
                    a2 = a_expand(index_a2);

                    % Calculate demand
%                     demandThisPeriod = demand(water, population(t), t, gwParam);
                    demandThisPeriod = gwParam.pumpingRate;

                    % Calculate cost and shortages this period
                    [ cost, ~,~, ~,~, ~, ~, ~, ~] = cost_supply_func( a1, a2, s1, s2, costParam, water, gwParam, t, demandThisPeriod, runParam.capacityDelay, exp_vectors, false);

                    % Calculate transition matrix
                    
                    % If stop pumping, move to state -1. Otherwise, use
                    % T_gw calculated above. 
                    
                    switch a1
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
                        if a2 == 1
                            T_exp_delay1_ind = 2;
                        elseif a2 == 2
                            T_exp_delay2_ind = 2;
                        end
                        temp_index = vectorIndex([T_exp_online_ind T_exp_delay1_ind T_exp_delay2_ind], {s_exp_on', s_exp_delay1', s_exp_delay2'});
                        exp_index = find(s_expand == temp_index);
                        T_expand(exp_index) = 1;
                    else
                        if a2 == 0
                            T_expand(index_s2) = 1; % Stay in current state
                        elseif a2 == 1
                            T_expand(index_s2 + 1) = 1; % Move up one state
                        elseif a2 == 2
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
                    
                    stateMsg = strcat('t=', num2str(t), ', s1=', num2str(s1), ', a1=', num2str(a1), ', s2=', num2str(s2), ', a2=', num2str(a2))
                    disp(stateMsg)
                   % Check if best decision
                    checkV = cost + expV;
                    if checkV < bestV
                        bestV = checkV;
                        bestX1 = a1;
                        bestX2 = a2;
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

    



end

