    %% calibrationMaster
% Brent Boehlert
% March 17, 2015

% This code calibrates in three ways: (1) daily 1990 to 1995 climate w/
% monthly 1990 to 1995 runoff (sub91 w/interpolated station data); (2)
% daily 1990 to 1995 climate w/ GRDC runoff (sub91 w/no station runoff
% data); and (3) montly 1990 to 1995 climate w/GRDC runoff (5 external
% Danube sub-basins)

clear
clc

cd('L:\02 Projects\WB Zimbabwe\src\CLIRUN calibration');

%% Step 1. Setup
setup

%% Step 2. Initialize matlabpools
switches.open = false;
if switches.open
  openPools
end

%% Step 3. Run Calibration
basins = 1:num.basins; 
data = double(data);

daysMat = repmat(days, num.basins, length(num.yrRange));

p = reshape(data(:, :, num.yrRange, 1), num.basins, []);
t = reshape(data(:, :, num.yrRange, 2), num.basins, []);
pet = reshape(data(:, :, num.yrRange, 3), num.basins, []).*daysMat;

calibStruct = calib_parallel(p, t, pet, runoffBC, basins);

filename = fullfile(home.outputs, 'Runoff', 'calibResults');
save(filename,'calibStruct');
disp( ['Saved ', filename] );

% Calibration performance to excel
precip = permute(data(:, :, num.yrRange, 1), [2 3 1]);
calib2Excel

%% Step 4. Graphics
qcGraphics

