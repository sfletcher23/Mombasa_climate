

%% Lat Lon
lat = [89:-2:1,-1:-2:-89];
lon = -179:2:179;
ind_lat_start = find(lat == -3);
ind_lat_end = find(lat == -5);
ind_lon_start = find(lon == 39);
ind_lon_end = find(lon == 41);

%% Temperature

load('CMIP5data.mat', 'tas_45')
tas_45 = tas_45(6:end);

t_45_mwache = single(zeros([length(tas_45) 2400]));
for i=1:length(tas_45)
    temp = tas_45{i}(ind_lat_start:ind_lat_end, ind_lon_start:ind_lon_end, :);
    t_45_mwache(i,:) = squeeze(mean(mean(temp,1),2))';
end
t_45_mwache = t_45_mwache([1:12,15:21,23:24],:)';
save('CMIP5_RCP45', 't_45_mwache')
clear tas_45 t_45_mwache

%% Precipitation

% load('CMIP5data.mat', 'pr_45')
% 
% pr_45 = pr_45(6:end);

p_45_mwache = single(zeros([length(pr_45) 2400]));
for i=1:length(pr_45)
    temp = pr_45{i}(ind_lat_start:ind_lat_end, ind_lon_start:ind_lon_end, :);
    p_45_mwache(i,:) = squeeze(mean(mean(temp,1),2))';
end
p_45_mwache = p_45_mwache([1:12,15:21,23:24],:)';
save('CMIP5_RCP45', 'p_45_mwache')
clear pr_45 p_45_mwache