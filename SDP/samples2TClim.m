
function [T_Precip, T_Temp] = samples2TClim(NUP, NUT, s_P, s_T, N)

P_step = s_P(2) - s_P(1);
T_step = s_T(2) - s_T(1);
P_bins = [s_P-P_step/2 s_P(end)+P_step/2];
T_bins = [s_T-T_step/2 s_T(end)+T_step/2];
M_P = length(s_P);
M_T = length(s_T);

% Check how many samples outside bins
numPSamp = numel(NUP);
numTSamp = numel(NUT);
indP = find(NUP < P_bins(1) | NUP > P_bins(end));
indT = find(NUT < T_bins(1) | NUT > T_bins(end));
percPout = length(indP)/numPSamp
percTout = length(indT)/numTSamp


% Calculate transition matrices
T_Precip = zeros(M_P, M_P, N);
T_Temp = zeros(M_T, M_T, N);
for t = 1:N-1
    for index_s_P = 1:M_P
        T_Precip(:,index_s_P,t) =  histcounts(NUP(:,t+1,index_s_P), P_bins, 'Normalization', 'Probability');
    end
    for index_s_T = 1:M_T
        T_Temp(:,index_s_T,t) =  histcounts(NUT(:,t+1,index_s_T), T_bins, 'Normalization', 'Probability');
    end
end



