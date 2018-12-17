function [runoff] = TP2runoff(T_ts, P_ts, steplen)

% Use hydrological model to esimate runoff from T and P timeseries

% Input
% T_ts is monthly T time series matrix in degrees C, numsamp x steplen*12
% P_ts is mohtly P time series matrix in mm/m, numsamp x steplen*12
% steplen is number of years in time series

% Output
% runoff is from CLIRUN in MCM/y , numsamp x steplen*12

% Load CLIRUN calibration results
calibrationFile = '29_Jan_2018_17_10_19_maybe_winner_3.mat';

if ~isempty(getenv('SLURM_JOB_ID'))
    load(strcat('/net/fs02/d2/sfletch/Mombasa_climate/CLIRUN/OutputData/data/',calibrationFile), 'X_results')
else
    load(strcat('/Users/sarahfletcher/Dropbox (MIT)/Mombasa_Climate/CLIRUN/OutputData/data/',calibrationFile), 'X_results')
end


numMonths = steplen*12;
[numSamples,~] = size(T_ts); 

streamflow_mmpd = zeros(numSamples,numMonths);
for i = 1:numSamples
    watyear = 1;
     streamflow_mmpd(i,:) = Simulator(X_results, T_ts(i,:), P_ts(i,:), watyear); % mm/d
end

area = 2250 * 1E6; % m2
streamflow_cmpd = streamflow_mmpd /1E3 * area;
streamflow_mcmpy = cmpd2mcmpy(streamflow_cmpd);
runoff = streamflow_mcmpy;

if false
    figure;
    plot(streamflow_mcmpy(1:10,:)')   
end


end