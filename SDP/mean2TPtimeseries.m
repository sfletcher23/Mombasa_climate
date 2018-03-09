function [T_ts, P_ts] = mean2TPtimeseries(timestep, steplen, sp, st, numsamp)

% Create timeseries of temperature and precipitation for a certain meam T
% and P states. 

% Inputs:
% timestep: refers to the time step of the SDP eg the 2nd two-decade period
% steplen: the length of the time step in years eg 20
% numsamp: number of time series samples to generate
% sp: mean precipitation value / state of SDP
% st: mean temperature value / state of SDP

% Outputs:
% T_ts: monthly temp time series (numsamp x timestep*12)
% P_ts: monthly precip time series (numsamp x timestep*12)


% Create index for appropriate dates in Tij and Pij
startdate = (1990-1900)*12 + 1; % First period starts in 1990
daterange = (timestep-1)*steplen*12 + startdate: (timestep)*steplen*12 + startdate - 1; % Second period starts in 2010

% Load Pij and Tij
if ~isempty(getenv('SLURM_JOB_ID'))
    load('/net/fs02/d2/sfletch/Mombasa_climate/BMA_code/Input/Mombasa_TandP.mat')
else
    load('/Users/sarahfletcher/Documents/MATLAB/Mombasa_Climate/BMA_code/Input/Mombasa_TandP.mat', 'Tij', 'Pij')
end

% Create standardized anomalies relative to decadal average across all GCMs for the current daterange
Tmean_gcm = mean(Tij(daterange,:),1);
Pmean_gcm = mean(Pij(daterange,:),1);
T_anoms = Tij(daterange,:) - repmat(Tmean_gcm, length(daterange),1);
P_anoms = (Pij(daterange,:) - repmat(Pmean_gcm, length(daterange),1)) ./ repmat(Pmean_gcm, length(daterange),1);

% Sample a GCM and one year of timeseries data, repeat x20 
% Do this numsamp times

model_rand = randi(21, numsamp, steplen);
year_rand = randi(steplen, numsamp, steplen);
time_ind_start = (year_rand - 1)*12 + 1;
time_ind_end = year_rand*12;

Tmean = st;
Pmean = sp;

T_ts = zeros(numsamp, steplen*12);
P_ts = zeros(numsamp, steplen*12);
for t = 1:steplen
    T_ts(:,(t-1)*12+1:t*12) = T_anoms(time_ind_start(:,t):time_ind_end(:,t) , model_rand(:,t))' + Tmean;
    P_ts(:,(t-1)*12+1:t*12) = P_anoms(time_ind_start(:,t):time_ind_end(:,t) , model_rand(:,t))' * Pmean + Pmean;
end

if false
    figure;
    subplot(2,1,1)
    plot(T_ts')
    subplot(2,1,2)
    plot(P_ts')
end




end