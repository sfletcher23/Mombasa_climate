%% Runoff time series for each T,P state

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
