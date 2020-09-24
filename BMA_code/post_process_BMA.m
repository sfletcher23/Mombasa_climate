%% Post processing BMA data

% This file takes the individual .csv files generated by the R code that
% does the Bayesian statistical analysis and integrates them into a single
% .mat file. 

% Inputs: .csv files of the form muUP_2020_scen1_job__2018-01-28_.csv

% Output: a .mat file with the following variables:

    % MUP: a numSamp x numTime x numStates (1000 x 5 x 31) matrix with
    % samples of the historical climate precipitation percent change
    
    % NUP: a numSamp x numTime x numStates (1000 x 5 x 31) matrix with
    % samples of the future climate precipitation percent change
    
    % MUT: a numSamp x numTime x numStates (1000 x 5 x 31) matrix with
    % samples of the historical climate temperature change
    
    % NUT: a numSamp x numTime x numStates (1000 x 5 x 31) matrix with
    % samples of the future climate temperature change
    
    % lamP: a numSamp x numTime x numStates (1000 x 5 x 31) matrix with
    % samples of precipitation reliabiltiies
    
    % lamT: a numSamp x numTime x numStates (1000 x 5 x 31) matrix with
    % samples of temperature reliabiltiies
    
%% Setup file reading  

% Which files to open
dateOpenT = '2019-07-12';
dateOpenP = '2019-07-09';
delta = false;
if delta
    outFolder = 'Output_deltas';
else
    outFolder = 'Output';
end

% Need to input jobid of T and P runs if running on cluster
jobid = (getenv('SLURM_JOB_ID'));
if ~isempty(jobid)
    jobIdT = num2str(67890);
    jobIdP = num2str(67889);
else
    jobIdT = '';
    jobIdP = '';
end

% How many P and T values

Nscens_T = 31;
Nscens_P = 31;

%% Loop over files and integrate into .mat file

% Temperature
for scen_ii = 1:Nscens_T
	scenind = scen_ii;
    
        time_ii = 0;
        for year = 2010:20:2090
            
            time_ii = time_ii+1;
            
            % Historical
            tmpstr = strcat(outFolder, '/', sprintf('muUT_%d_scen%d',year,scen_ii),'_','job', '_',jobIdT,'_', dateOpenT,'_.csv');
            tmp = csvread(tmpstr);
            MUT(:,time_ii,scenind) = tmp;
            
            % Future
            tmpstr = strcat(outFolder, '/',sprintf('nuUT_%d_scen%d',year,scen_ii),'_','job', '_',jobIdT,'_', dateOpenT,'_.csv');
            tmp = csvread(tmpstr); 
            NUT(:,time_ii,scenind) = tmp; 
            
            % Lambda
            tmpstr = strcat(outFolder, '/',sprintf('lambdaUT_%d_scen%d',year,scen_ii),'_','job', '_',jobIdT,'_', dateOpenT,'_.csv');
            tmp = csvread(tmpstr);
            lamT(:,:,time_ii,scenind) = tmp;
	     
        end
end

% Precipitaiton
for scen_ii = 1:Nscens_P
    scenid = scen_ii;

    time_ii = 0;
    for year = 2010:20:2090
        
        time_ii = time_ii+1;

        % Historical
        tmpstr = strcat(outFolder, '/',sprintf('muUP_%d_scen%d',year,scen_ii),'_','job', '_',jobIdP,'_', dateOpenP,'_.csv');
        tmp = csvread(tmpstr);
        MUP(:,time_ii,scenid) = tmp;

        % Future
        tmpstr = strcat(outFolder, '/',sprintf('nuUP_%d_scen%d',year,scen_ii),'_','job', '_',jobIdP,'_', dateOpenP,'_.csv');
        tmp = csvread(tmpstr);
        NUP(:,time_ii,scenid) = tmp;

        % Lambda
        tmpstr = strcat(outFolder, '/',sprintf('lambdaUP_%d_scen%d',year,scen_ii),'_','job', '_',jobIdP,'_', dateOpenP,'_.csv');
        tmp = csvread(tmpstr);
        lamP(:,:,time_ii,scenid) = tmp;
    end
end    

% Save
saveName = strcat('BMA_results_', datestr(datetime, 'yyyy-mm-dd')); 
save(saveName,'MUP','MUT','NUT','NUP','lamT','lamP', 'jobIdT', 'jobIdP', 'delta')



