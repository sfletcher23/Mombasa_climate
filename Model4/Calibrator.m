addpath('Functions');
calibYears = {'19716' '2000'}; 
MTHavgCalib = 0;
FN = 'UNH_All22';
dataFN = 'inputdata/Vietnam_Hydro_Observed.mat';
init = [  100   50  0.5   0.001   .1    1.1  .2  ]; % initial conditions
init = [  10.57   157.98  0.7637   0.001   .101    1.099  .2  ]
%init =[11.5694964599609,151.083862762451,0.76760625,0.001,0.0384904632568359,1.09484960937500,0.23125];
filename = {[FN,'_1']};
CLIRUN_II_calibrator_obs(filename,dataFN,init,calibYears,MTHavgCalib)
clear
