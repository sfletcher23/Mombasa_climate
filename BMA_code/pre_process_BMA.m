
% BMA pre-processing

% This file prepares historical climate data and climate model projections
% for input into the Bayesian statistical model

% Input data is available in the file Input/Mombasa_TandP.mat and includes
% the following:

    % Historical precipitation (P0) and temperature (T0) data is sourced from CRU dataset version TS.3.21
        % Monthly data from 1901 to 2012 
        % Averaged over 2 degrees S to 6 degrees S and 38 degrees E to 42 degrees E

    % An ensemble of 21 GCM projections of precipitaiton (Pij) and temperature (Tij)
        % Models listed in SI Table 1, regridded as described in methods.
        % Monthly values from 1901-2100
        
% Saved output data in .csv files includes:

    % X_year files, where year ranges from 1990 to 2090
        % Contains 20-year average from 10 years before to 10 years after.
        % Eg. 1990 file contains data averaged from 1980 to 2000.
        % Two rows: top row is T, bottom row is log P
        % 21 columns corresponding to climate model ensembles
        
    % Pij and Tij saved as .csv files with 2400 rows for each monthly value
    % and 21 columns for the climate ensembles
    
    % X0PU and X0TU are .csv files that contain the virtual future
    % observations of climate that will be used to condition ???
        


%% Load data and setup 

% Load data
clear all; close all
load('Input/Mombasa_TandP.mat'); 

% Process historical and projected T and P data and save in .csv files
rewrite_XYlambda = true;

% Create virutal future observations of T and P to condition on
create_scens = true;


%% Creating initial X, Y and lambda values
if rewrite_XYlambda
    
    % Select historical data from 1970 to 2000
    % Transform P to log
    for year = 70:100
        YT0(year-69) = mean(T0(12*(year-1)+1:12*year));
        YP0(year-69) = log(mean(P0(12*(year-1)+1:12*year)));
    end
    
    % Use historical mean and std to develop prior dist for lambda
    % (reliabilities)
    X0 = [mean(YT0), mean(YP0)]';
    SD = [std(YT0), std(YP0)]';
    lambda0 = SD.^(-2);
    
    % Take yearly averages for projection data
    % Transform P to log
    for year = 1:200
        YTij(year,:) = mean(Tij(12*(year-1)+1:12*year,:),1);
        YPij(year,:) = log(mean(Pij(12*(year-1)+1:12*year,:),1));
    end

    % Saving variables for R code: 
    csvwrite('Input/lambda0.csv',lambda0)
    csvwrite('Input/Tij.csv',Tij)
    csvwrite('Input/Pij.csv',Pij)
    
    % Take decadal averages from 1990 to 2090
    for decade = 1:11
        X = [mean(YTij(10*(decade-1)+80:10*(decade-1)+100,:),1)', mean(YPij(10*(decade-1)+80:10*(decade-1)+100,:),1)']';
        str1 = sprintf('Input/X_%2.0f.csv',1990+10*(decade-1));
        csvwrite(str1,X)
    end
    
end

%% Creating virtual future observations of T and P

% Here we create the range and discretiztation of future observations of T
% and P that are used to condition the Bayesian model. 

% These may not be the actual values, but used for exploratory analysis ???

if create_scens
    deltaT_X0(1:10,1) = zeros(10,1);
    deltaT_X0(1:10,2) = 0:0.25:0.25*9;
    deltaT_X0(1:10,3) = 0:0.5:4.5;

    deltaP_X0(1:10,1) = 0:-0.08:-0.72;
    deltaP_X0(1:10,2) = zeros(10,1);
    deltaP_X0(1:10,3) = 0:0.08:0.72;

    % Univariate case
    X0_Topts = X0(1) + deltaT_X0;
    X0_Popts = log(exp(X0(2))+exp(X0(2))*deltaP_X0);

    csvwrite('Input/X0TU.csv',X0_Topts)
    csvwrite('Input/X0PU.csv',X0_Popts) 
end
