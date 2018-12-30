

% BMA pre-processing

% This file prepares historical climate data and climate model projections
% for input into the Bayesian statistical model

% Input data is available in the file Input/Mombasa_TandP.mat and includes
% the following:

    % Historical precipitation (P0) and temperature (T0) data is sourced from CRU dataset version TS.3.21
        % Monthly data from 1901 to 2012 
        % Averaged over 2 degrees S to 6 degrees S and 38 degrees E to 42 degrees E
        % Precipitation is in mm/month, temperature is in degrees C

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
    % observations of climate that will be used to condition the
    % uncertainty estimates in the next time period
        


%% Creating initial X, Y and lambda values.  These don't change 

load('Mombasa_TandP.mat'); 


% I think this commented section below is not used in R code. Okay to delete?
% % Select historical data from 1960 to 2010
% % Log transform on P data so closer to normally distributed
% for year = 61:110
%     YT0(year-60) = mean(T0(12*(year-1)+1:12*year));
%     YP0(year-60) = log(mean(P0(12*(year-1)+1:12*year)));
% end

% % Calculate deltas: change between means in 20-year consecutive time
% % periods
% X0 = [mean(YT0(21:40))-mean(YT0(1:20)), mean(YP0(21:40))-mean(YP0(1:20))]';

% Take yearly averages from projection data
for year = 1:200
    YTij(year,:) = mean(Tij(12*(year-1)+1:12*year,:),1);
    YPij(year,:) = mean(Pij(12*(year-1)+1:12*year,:),1);
end

% I think this commented section below is not used in R code. Okay to delete?
% % Saving variables for R code:   
% csvwrite('Input/Tij.csv',Tij)
% csvwrite('Input/Pij.csv',Pij)

% Loop over 6 time periods (20 years each, first one is 1960-1980)
for dec = 1:6
    
    % Calculate change in average temperature from one period to next
    X(1,:) = [mean(YTij(20*(dec-1)+81:20*(dec-1)+100,:),1)'-mean(YTij(20*(dec-1)+61:20*(dec-1)+80,:),1)']';
    
    % Calculate change in average precipitation from one period to next
    X(2,:) = [(mean(YPij(20*(dec-1)+81:20*(dec-1)+100,:),1)'-mean(YPij(20*(dec-1)+61:20*(dec-1)+80,:),1)')./mean(YPij(20*(dec-1)+61:20*(dec-1)+80,:),1)']';
    
    % Save data for historical time period
    str1 = sprintf('Input/X_%2.0f.csv',1990+20*(dec-1));
    csvwrite(str1,X)
    
    % Calculate lambda priors using std dev of climate data
    SD = [std(X(1,:)), std(X(2,:))]';
    lambda0 = SD.^(-2);
    str1 = sprintf('Input/lambda0_%2.0f.csv',1990+20*(dec-1));
    csvwrite(str1,lambda0)
end

%% Creating files with virtual future climate observations

deltaT_X0 = repmat([0:0.05:1.5],10,1);
deltaP_X0 = repmat([-0.3:0.02:0.3],10,1);

% Univariate case
X0_Topts = deltaT_X0;
X0_Popts = deltaP_X0;

csvwrite('Input/X0TU.csv',X0_Topts)
csvwrite('Input/X0PU.csv',X0_Popts) 

