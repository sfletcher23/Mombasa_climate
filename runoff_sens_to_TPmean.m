% Test sensititivty of mean T and P on runoff 

% Need to run first few sections of climate sdp first for T,P time series
% and parameters


T_ts_base = T_ts{120,5};
P_ts_base = P_ts{30,5};
runoff_base = TP2runoff(T_ts_base, P_ts_base, runParam.steplen);
MAR_base = mean(runoff_base)

T_ts_sens = T_ts_base - climParam.T_delta/2;
P_ts_sens = P_ts_base + climParam.P_delta/2;
runoff_sens = TP2runoff(T_ts_sens, P_ts_sens, runParam.steplen);
MAR_sens = mean(runoff_sens)

MAR_diff = MAR_sens - MAR_base
MAR_perc = MAR_diff / MAR_base