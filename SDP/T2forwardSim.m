
function [T_over_time, P_over_time] = T2forwardSim(T_Temp, T_Precip, s_P, s_T, N, t_now, T0, P0, numSamp)

p = rand(numSamp,N);
state_ind_P = zeros(numSamp,N-t_now);
state_ind_T = zeros(numSamp,N-t_now);
state_ind_P(:,1) = find(P0==s_P);
state_ind_T(:,1) = find(T0==s_T);
for i = 1:numSamp
    for t = 1:N-t_now
        time = t_now + t -1;
        state_ind_T(i,t+1) = find(p(i,t) < cumsum(T_Temp(:,state_ind_T(i,t),time)),1);
        state_ind_P(i,t+1) = find(p(i,t) < cumsum(T_Precip(:,state_ind_P(i,t),time)),1);
    end
end
T_over_time = s_T(state_ind_T);
P_over_time = s_P(state_ind_P);

