function [T_ts, P_ts] = mean2TPtimeseries(T, P, NUP, NUT, timestep, steplen, index_s_p, index_s_t, numsamp)

% Create timeseries of temperature and precipitation for a certain meam T
% and P states. 

% Inputs:
% timestep: refers to the time step of the SDP eg the 2nd two-decade period
% steplen: the length of the time step in years eg 20
% index_s_p: the index of the precipitation state, which corresponds to the third dimension of NUP
% index_s_t: the index of the temperature state, which corresponds to the third dimension of NUT
% NUT: BMA mean temp samples (samples x timesteps x numstates)
% NUP: BMA mean precip samples (samples x timesteps x numstates)

% Create index for appropriate dates in Tij and Pij
startdate = (1990-1900)*12 + 1; % First period starts in 1990
daterange = (timestep-1)*steplen*12 + startdate: (timestep)*steplen*12 + startdate - 1; % Second period starts in 2010

% Load Pij and Tij
load('/Users/sarahfletcher/Documents/MATLAB/Mombasa_Climate/BMA_code/Input/Mombasa_TandP.mat', 'Tij', 'Pij')

% Create standardized anomalies relative to decadal average across all GCMs for the current daterange
Tmean_gcm = mean(Tij(daterange,:),1);
Pmean_gcm = mean(Pij(daterange,:),1);
T_anoms = Tij(daterange,:) - repmat(Tmean_gcm, length(daterange),1);
P_anoms = (Pij(daterange,:) - repmat(Pmean_gcm, length(daterange),1)) ./ repmat(Pmean_gcm, length(daterange),1);

% Sample a GCM and one year of timeseries data, repeat x20 
% Do this numsamp times

[bma_samp_size,~,~] = size(NUT);
samp_rand = randi(bma_samp_size, 1, numsamp);
model_rand = randi(21, numsamp, steplen);
year_rand = randi(steplen, numsamp, steplen);
time_ind_start = (year_rand - 1)*12 + 1;
time_ind_end = year_rand*12;

Tmean = NUT(samp_rand, timestep, index_s_t);
Pmean = NUT(samp_rand, timestep, index_s_p);

T_ts = zeros(numsamp, steplen*12);
P_ts = zeros(numsamp, steplen*12);
for t = 1:steplen
    T_ts(:,(t-1)*12+1:t*12) = T_anoms(time_ind_start(:,t):time_ind_end(:,t) , model_rand(:,t))' + repmat(Tmean,1,12);
    P_ts(:,(t-1)*12+1:t*12) = P_anoms(time_ind_start(:,t):time_ind_end(:,t) , model_rand(:,t))' .* repmat(Pmean,1,12) + repmat(Pmean,1,12);
end





end