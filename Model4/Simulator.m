%CLIRUN SIMULATOR
clear all;
global  PET NYRS WATYEAR % DATA INPUT
global  PRECIP_0% CALCULTED VALUES daily
global  precip_day days % CALCULTED VALUES 
global ku  kp  kl lm  sat over

addpath('Functions')

senar={'bccr_bcm2_0','cccma_cgcm3_1','cccma_cgcm3_1_t63','cnrm_cm3' ...
'csiro_mk3_0','csiro_mk3_5','gfdl_cm2_0','gfdl_cm2_1','giss_aom' ...
'giss_model_e_h','giss_model_e_r','iap_fgoals1_0_g','inmcm3_0','ipsl_cm4','miroc3_2_hires' ...
'miroc3_2_medres','mpi_echam5','mri_cgcm2_3_2a','ncar_ccsm3_0' ...
'ncar_pcm1','ukmo_hadcm3','ukmo_hadgem1'}

dispOn = 1;
makePet = 1;

%deleet 
senar={'Vietname_Hydro_UNH'}

%

for(ind=1:1)

 FN1 = senar{ind};
for(Senario_Number=1:1)
%%--------------------------------------------
%days in each month
        WATYEAR = 6;
        days = [31 28 31 30 31 30 31 31 30 31 30 31];% CRU -JAN- DEC ;
%         WATYEAR = 10;
%         days = [31 30 31 31 28 31 30 31 30 31 31 30 ];% Oct- Sept;

% load data
load(['InputData/ForSims/',FN1,'S',num2str(Senario_Number,1),'.mat'], '-regexp',... 
    'PRECIPTS', 'TEMPTS', 'TRANGETS', 'LATm', 'flowNames') % historic weather data and obs flows

start = 1901;

load('InputData/16_Jan_2018_14_20_14_Pton_point_normparam.mat', '-regexp', 'X_results')
        
nbasins = size(TEMPTS,1);
basins = 1:nbasins; %*
%basin2clim = (1:nbasins); % for CRU - basins and station match


%%
%set up simulation time limits

sessionName = FN1; 
answer = {'1901';'2000'};
BASIN = ['OutputData/data/sim/',sessionName,'_','S',num2str(Senario_Number,1),'_',answer{1},'to',answer{2}];
display(BASIN);


d = str2double(answer);
startSim = (d(1));
finishSim = (d(2));


NYRS = finishSim - startSim +1 ;
months = NYRS*12;
timeOne = (startSim-start)*12+1;
timeFrame = timeOne:(timeOne+months-1);

%%
if makePet
    %calculate PET for all CRU locations
    pet= zeros(size(TRANGETS));
    for loc = 1:nbasins
        
        %change units
        %PRECIPTS = double( precip(loc,:) )/10;
        %TEMPTS = double( temp(loc,:) )/10;
        %TRANGETS = double( trange(loc,:) )/10;
        
        LAT = LATm(loc);
       
        TEMP = TEMPTS(loc,:);
        TRANGE = TRANGETS(loc,:);
        PRECIP = PRECIPTS(loc,:);
        
        %Call PET model
        PET= ModHargreaves3(LAT,WATYEAR,TEMP,TRANGE,PRECIP);
        
        pet(loc,:) = PET;
    end
end

%% Begin simulation

ROresults=zeros(nbasins,months);

for bas=1:nbasins
    %     
    PET = [ pet(bas,timeFrame) , 0 ];
    PRECIP_0 = [ PRECIPTS(bas,timeFrame)./repmat(days,1,NYRS) , 0 ];% make mm/day
    
    %Message
    disp('SIMULATING...');
    disp(['> Basin: ',num2str(basins(bas)),'  ',flowNames{bas}]);
    disp(['> From: ',num2str(startSim)]);
    disp(['> To: ',num2str(finishSim)]);
    disp(['> Started At: ',datestr( now )]);

    % bring in calibration results
    
    x = X_results(bas,:);
    
    %set calib params for ode45 of 'flowul01nile'
    %layer thickness
    sat= x(1);  lm = x(2);
    % flux parameters
    ku= x(3);kp = x(4);kl= x(5);
    % coef
    inter = x(6); over = x(7);

    precip_day = inter*PRECIP_0;
    
    Tspan= (0:months)';
    options = odeset('NonNegative',[1 2]);
    [Tvec,xsol]=ode45('soil_model_II',Tspan,[5,.1,0.0,0.0], options);%JMS-mod**
    
    
    % Conforming raw data from ODE
    z = xsol(2:end,1:2);
    wb = runoff_II(z);
    wb = wb';
    
%if dispOn
 %   figure
 %   plot(wb(:,4), 'DisplayName', 'xout(1:180,1:2)', 'YDataSource', 'xout(1:180,1:2)'); figure(gcf)
%end
    
ROresults(bas,:) = squeeze(wb(:,4)');
end

save([BASIN,'.mat'],'ROresults','basins','startSim','finishSim')
%%--------------------------------------------
end 
end
% RO=ROresults;
% load M_hist_Sim_1961_1990 ROresults;
% delRO=RO-ROresults; ann_delRO=sum(delRO,2)/30;
% ann_histRO=sum(ROresults,2)/30; Nann_delRO=100*ann_delRO./ann_histRO;
% A=[sub_basin ann_delRO Nann_delRO];
% csvwrite('data/M_GFDLcm21_delROann.csv', A);

