% setup 

%% Step 1. Home
home.home = 'L:';
home.data = [home.home '\01 Data'];
home.proj = [home.home, '\02 Projects\WB Zimbabwe'];
home.inputs = [home.proj, '\Inputs'];
home.outputs = [home.proj, '\Outputs'];
home.clirun = [home.home, '\01 Data\CLIRUN\clirun'];

% Paths
addpath([home.proj, '\src']);
addpath(genpath([home.home, '\01 Data\src']));
addpath([home.clirun, '\model']);
addpath([home.clirun, '\snow']);

%% Step 2. Load calibration data
load([home.inputs '\MAT files\clirunCalInputs']);  
num.yrRange = 14:33;  %1961 to 1980
num.months = 12;

names.basins = key.basins;
num.basins = length(names.basins);

days = [31 28 31 30 31 30 31 31 30 31 30 31];

clear key
