
% CALIBRATION RESULTS

postProcess = false;
if postProcess 
    
    global PRECIP_0 OBS NYRS PET % DATA INPUT
    global precip_day % CALCULTED VALUES
    global ku kp kl sat lm inter over % MODEL PARAMETERS
    
    usePton = false;
    MTHavgCalib = 0; % 1 monthly mean, 0 individual points
    
    load('InputData/Mwache_hydro_data')
    x = X_results;
    bas = 1;
    
    days = [31 28 31 30 31 30 31 31 30 31 30 31];% CRU -JAN- DEC ;
    WATYEAR = 6;   % Which year of the month start on
    days_adj = [days(WATYEAR:end) days(1:WATYEAR-1)];
    NYRS = 14;
    months = NYRS*12;
    makePet = false;
    PET = PET_date_lin'; % mm/day
    PRECIP_0 = P0_date_lin' ./ repmat(days_adj,1,NYRS); % Convert from mm/m to mm/day
    OBS = runoff_mwache; % mm/day
    flowNames = {'Mwache river'};


    if usePton
       load('InputData/Mwache_hydro_pton') 
       PET = E_pton; % mm/day
       PET(118) = (PET(119) + PET(120)) /2;
       PRECIP_0 = P_pton; % mm/day
       OBS = runoff_pton; % mm/day
       OBS(118) = (OBS(119) + OBS(120)) /2;
       PRECIP_0(118) = (PRECIP_0(119) + PRECIP_0(120)) /2;
    end
end


%set globals for ode45 of 'flowul01nile'
%model parameters from OPTIMIZATION CODER
%layer thickness
sat= x(1);
lm = x(2);
% flux parameters
ku= x(3);
kp = x(4);
kl= x(5);
% coef
inter = x(6);
over = x(7);
%set precip data
    
precip_day = inter.*PRECIP_0;

months = NYRS*12;

Tspan= (0:months)';

options = odeset('NonNegative',[1 2]);
[Tvec,xsol]=ode45('soil_model_II',Tspan,[5,.1,0.0,0.0], options);%JMS-mod
    
    % Conforming raw data from ODE
    z = xsol(2:end,1:2);
    wb = runoff_II(z);
    wb = wb';
    
    if MTHavgCalib==2
        xout = zeros(12,2);
        xout(:,1) = mean( reshape( wb(:,4)', 12, []), 2)';
    else xout = zeros(length(wb),2);
    % ODE SOVLERS OUTPUT TIME O - STRIPPING OFF TIME ZERO 1 to 12 * years
    % wb = [rss1; dr1; rs1; RunOff];
    xout(:,1) = wb(:,4)'; % RunOff - total runoff
    end
    
    xout(:,2) = OBS;
    

    % ESTABLISH THE OBJECTIVE FUNCTION
    
    %To find month difference model v. observed
    ken  = xout(:,2) - xout(:,1);
    diff2 = ken.^2;
    
    tssum   = sum(diff2);
    obj =  tssum;
    
    if isreal(obj);
    
%     display(obj)
    
    XX = sum(xout(:,1)); % Model
    YY = sum(xout(:,2)); % Obs
    error = (YY-XX)/YY;  % Obs - Model/
    r2 = 1 - (tssum ./ sum((xout(:,2) - mean(xout(:,2))).^2));
    disp('~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-')
    disp(['r2 (point): ',num2str(r2)])
    disp(['Error: ',num2str(error)])
    disp('~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-')
    
        %to plot..TOTAL runoff v. Observes         %JMS
   figure;
   subplot(2,1,1)
   plot(xout(:,1:2), 'LineWidth', 1),  legend({'Mdl R','obs R'}, 'Location', 'northwest'); %JMS
   xlabel('month'), ylabel('runoff')
   legend('boxoff')
   ylim([0 3.5])
   
   % Do same error analysis with monthly means instead of individual points
   xout1_mean = mean( reshape( xout(:,1),12,[]) , 2 )';%Make monthly mean Model
   xout2_mean = mean( reshape( xout(:,2),12,[]) , 2 )';% 12 monthlty mean obs
   
   ken  = xout2_mean - xout1_mean;
    diff2 = ken.^2;
    tssum   = sum(diff2);
    obj =  tssum;
    r2 = 1 - (tssum ./ sum(( xout2_mean - mean( xout2_mean)).^2));
    
    XX = sum(xout1_mean); % Model
    YY = sum(xout2_mean); % Obs
    error = (YY-XX)/YY;  % Obs - Model/
    r2 = 1 - (tssum ./ sum((xout2_mean - mean(xout2_mean)).^2));
    disp('~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-')
    disp(['r2 (mean): ',num2str(r2)])
    disp(['Error: ',num2str(error)])
    disp('~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-')
    
    subplot(2,1,2)
   plot(1:12, [xout1_mean; xout2_mean], 'LineWidth', 1),  legend({'Mdl mean R','obs mean R'},'Location', 'northwest'); %JMS
   xlabel('month'), ylabel('runoff')
   legend('boxoff')
   ylim([0 0.8])
   
   ans = X_results'

    else display(bas); display('obj not real');
    end
    
    if ~postProcess
    saveas(1,['OutputData/plots/calib_','_basin',num2str(bas),'_',flowNames{bas},'.jpg'])
    
    RESULTS = ones(1,3);
    RESULTS(1) = obj;
    RESULTS(2) = error;  
    RESULTS(3) = r2;
    
    
    
    save(['OutputData/data/Calibration_','_',flowNames{bas},'.mat'],'RESULTS','x','xout')

    X_results(bas,:)=x;
    resultsMAT(bas,:)=RESULTS;
    %xoutMAT(bas,1:months,:)=xout;
    xoutMAT(bas,:,:)=xout;
    end
    
  