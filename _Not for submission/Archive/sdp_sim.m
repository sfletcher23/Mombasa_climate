function [ C_state, T_state, P_state, action, damCostTime, shortageCostTime,...
    opexCostTime, totalCostTime] = sdp_sim( climParam, runParam, costParam, s_T_abs, s_P_abs, T_Temp, T_Precip, s_C, M_C, X, shortageCost, a_exp, infra_cost, desal_opex)

R = 1000;
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


for i = 1:R
    for t = 1:N
        
        % Choose best action given current state
            index_t = find(T_state(i,t) == s_T_abs);
            index_p = find(P_state(i,t) == s_P_abs);
            
        
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
        
        for k = 1:4
            
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
            shortageCostTime(i,t,k) = shortageCost(index_t, index_p, short_ind, 1)  / (1+costParam.discountrate)^((t-1)*runParam.steplen+1);
            if t == 1 
                shortageCostTime(i,t,k) = 0;  % This is upfront building period
            end
            ind_dam = find(a == a_exp);
            damCostTime(i,t,k) = infra_cost(ind_dam)  / (1+costParam.discountrate)^((t-1)*runParam.steplen+1);
            if runParam.desalOn
                opexCostTime(i,t,k) = desal_opex(index_t, index_p, short_ind, t);
            end
            totalCostTime(i,t,k) = (shortageCostTime(i,t,k) + damCostTime(i,t,k));
            totalCostTime(i,t,k) = totalCostTime(i,t,k) + opexCostTime(i,t,k);
            
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