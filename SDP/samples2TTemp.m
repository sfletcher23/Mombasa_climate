
function [ T_Temp] = samples2TTemp( NUT, s_T, N)

T_step = s_T(2) - s_T(1);
T_bins = [s_T-T_step/2 s_T(end)+T_step/2];
M_T = length(s_T);

% Check how many samples outside bins

numTSamp = numel(NUT);
indT = find(NUT < T_bins(1) | NUT > T_bins(end));
percTout = length(indT)/numTSamp

% Calculate transition matrices

T_Temp = zeros(M_T, M_T, N);
for t = 1:N
    for index_s_T = 1:M_T
        T_Temp(:,index_s_T,t) =  histcounts(NUT(:,t,index_s_T), T_bins, 'Normalization', 'Probability');
    end
end



