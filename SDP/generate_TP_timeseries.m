% Generate T and P timeseries based on T and P mean states


N = 5;
climParam.T_delta = 0.05; % deg C
climParam.T0_abs = 26;
T_abs_max = max(s_T) * N;
s_T_abs = climParam.T0_abs : climParam.T_delta : climParam.T0_abs+ T_abs_max;
climParam.numSampTS = 1;
s_P_abs = 66:1:97;

T_ts = cell(M_T_abs,N);
P_ts = cell(M_P_abs,N);

for t = 1:N
    
    for i = 1:M_T_abs  
        [T_ts{i,t}, ~] = mean2TPtimeseriesMJL(t, runParam.steplen, s_P_abs(1), s_T_abs(i), climParam.numSampTS);
    end
    
    for i = 1:M_P_abs  
        [~, P_ts{i,t}] = mean2TPtimeseriesMJL(t, runParam.steplen, s_P_abs(i), s_T_abs(1), climParam.numSampTS);
    end
    
end
