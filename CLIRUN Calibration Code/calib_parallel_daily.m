function calibStruct =  calib_parallel_daily(PRECIP_TS,TEMP_TS,PET_TS,runoff,basins)
%function calibStruct =  calib_parallel_daily(PRECIP_TS,TEMP_TS,PET_TS,runoff,basins)
%
% Code assumes that pet, t, and p are at a daily time step, and runoff is
% at a monthly time step.  Modeled runoff is daily, but is summed to
% monthly to evaluate the objective function, which is minimizing the sum
% of squared errors between monthly observed and modeled runoff.

%NOTE ABOUT DAILY VERSUS MONTHLY
% Runoff inputs are assumed to be monthly, and are
% maintained as such for use in the objective function. 
% Precip and PET are assumed to be daily.

%days in each month
days = [31 28 31 30 31 30 31 31 30 31 30 31];

%set calibration variables
%model parameters from OPTIMIZATION CODER
%tl= x(1);th= x(2); ku= x(3);kp =x(4);kl=x(5);sat =x(6);
%inter = x(7); over=x(8); dm =x(9);
%init = [0 5 0.5  0.001 0.1 10 0.9 0.1 10];
init = [1 2 0.5  0.001 0.1 10 1.5 0.1 10];
lm=100;

%set upper and lower bounds
%up=    [ 25  50 0.9  .1      .5   90  1.0   0.3   50];
%lo=    [-25   0 0.01 .0001  0.001  2   .8   0.01   0];
up=    [ 5  10 0.9  .1      .5   90  3.0   0.3   50];
lo=    [-10  0 0.01 .0001  0.001  2   .8   0.01   0];

%linear constraints
AMat = [1 -1 0 0 0 0 0 0 0];
bVec = 0;

%%
[~, numDays]= size(PRECIP_TS);
[numBasins numMonths]= size(runoff);
numYears = numDays/365; %numMonths/12;
numMonthsRun = size(runoff, 2);
numYearsRun = numMonthsRun/12;

%Begin calibration
startTime = now;
calibStruct = repmat(struct('x',zeros(1,9),...
    'tl',0,'th',0,'ku',0,'kp',0,'kl',0,'sat',0,'inter',0,'over',0,'dm',0,...
    'xlsRow','',...
    'xlsResultsRow','',...
    'model',zeros(numDays, 1),...
    'modelMo', zeros(numMonths, 1),...
    'observed',zeros(numMonthsRun+2,1),...
    'obj',0,'slope',0,'r2',0,'nash',0,'error',0),...
    numBasins,1);
parfor bas=1:numBasins
    basinStr = num2str( basins(bas) );
    
    OBS = runoff(bas,:); %./repmat(days, 1, numYearsRun); %convert runoff to daily
    if any(isnan(OBS))
        display('CALIBRATION SKIPPED');
        display(['> Basin: ',basinStr]);
        display(' ');
        continue %skip if missing Data
    end
    PET = double([ PET_TS(bas,:) , 0 ]); %./repmat(days, 1, numYears) , 0 ]; %PET to daily
    TEMP = double([TEMP_TS(bas,:) , 0 ]);
    PRECIP_0 = double([ PRECIP_TS(bas,:) , 0 ]); %./repmat(days, 1, numYears) , 0 ]; %Precip to daily
    
    %options = psoptimset('OutputFcn', @myfun);
    options = psoptimset('CompletePoll','on', 'MaxIter', 100,'Display','off');
    
    %Message
    disp('CALIBRATING...');
    disp(['> Basin: ',basinStr]);
    disp(['> Started At: ',datestr( now ),]);
    disp(' ');
    %Run Pattern Search
    tic;
    twoLayerHandle = @(x)twolayer01nile_parallel_daily(x,PRECIP_0,OBS,TEMP,PET,numYears,lm);
    [x,fval,exitflag,output] = patternsearch(twoLayerHandle,...
        init, AMat, bVec, [], [], lo, up,[], options);
    %[x,fval,exitflag,output] = ga(@twolayer01nile_test2,...
    %    length(init), AMat, bVec, [], [], lo, up );
    
    calibStruct(bas).xlsRow = num2str(bas+1);
    
    display('CALIBRATION FINISHED');
    disp(['> Finished Basin ',basinStr,' at ',datestr( now )]);
    disp(['> Time Elapsed: ',num2str(toc),'sec']);
    totalTime = datevec(now - startTime);
    disp(['> Total Time Elapsed: ',num2str(totalTime(3)),':',...
        num2str(totalTime(4)),':',num2str(totalTime(5)),':',num2str(totalTime(6))]);
    disp(' ');
    %estCompletionTime = startTime + numBasins*(now - startTime)/(bas-basinsStart+1);
    %disp(['> Estimated Completion Time: ',datestr(estCompletionTime)]);
    
    
    % CALIBRATION RESULTS
    
    %set values for ode45 of 'flowul01nile'
    calibStruct(bas).x = x;
    tl =    x(1);
    th =    x(2);
    ku=     x(3);
    kp =    x(4);
    kl =    x(5);
    sat =   x(6);
    inter = x(7);
    over =  x(8);
    dm =    x(9);
    calibStruct(bas).tl =    x(1);
    calibStruct(bas).th =    x(2);
    calibStruct(bas).ku=     x(3);
    calibStruct(bas).kp =    x(4);
    calibStruct(bas).kl =    x(5);
    calibStruct(bas).sat =   x(6);
    calibStruct(bas).inter = x(7);
    calibStruct(bas).over =  x(8);
    calibStruct(bas).dm =    x(9);
    
    %set precip data
    %interMo = [inter inter inter ones(1,6) inter inter inter];
    PRECIP_DAY = inter*PRECIP_0;
    %PRECIP_MO = [repmat(interMo,1,numYears) 1].*PRECIP_0.*[repmat(days,1,numYears) 0];
    snowM = snowmodel_New(TEMP,PRECIP_DAY,tl,th,dm);
    precipStruct = snowMeltData(PRECIP_DAY, [0,snowM]);
    PRECIP = precipStruct.precipAvail;
    ISMELT = (PRECIP>PRECIP_DAY);
    precip_day = PRECIP; %./[repmat(days,1,numYears) 10];
    
    
    Tspan= (0:numDays)';
    options = odeset('NonNegative',[1 2]);
%     flowHandler = @(t,x)flowul01nile_parallelx(t,x,precip_day,PET,ku,kp,kl,lm,sat,over);
    flowHandler = @(t,x)flowul01nile_parallel(t,x,precip_day,PET,ku,kp,kl,lm,sat,over);
    [Tvec,xsol]=ode23(flowHandler,Tspan,[5,.1,0.0,0.0],options);%JMS-mod
    
    % Conforming raw data from ODE
    z = xsol(2:end,1:2);
    wb = runoff01nile_parallel_daily(z,numYears,precip_day,ku,kl,sat,over);
    wb = wb';
    
    
    % ODE SOVLERS OUTPUT TIME O - STRIPPING OFF TIME ZERO 1 to 12 * years
    % wb= [rss1; dr1; rs1; RunOff];
    
    %if numMonthsRun == 12;
    %  calibStruct(bas).model = mean( reshape( wb(:,4),12,[]) , 2 ); % RunOff - total runoff
    %else
      calibStruct(bas).model = wb(:,4); %BBB model to accommodate more months
    %end
    
    tempModel = reshape(calibStruct(bas).model, 365, []);
    tempModelMo = zeros(12, numYears);
    for mo = 1:12
      tempRange = 1+sum(days(1:mo))-days(mo):sum(days(1:mo));
      tempModelMo(mo, :) = sum(tempModel(tempRange, :), 1);
    end
    
    %if obs runoff is GRDC (or just one year), average modeled to 12 values
    if length(OBS) <= 12   
      tempModelMoRes = mean(tempModelMo, 2);
    else
      tempModelMoRes = reshape(tempModelMo, [], 1);
    end
    
    calibStruct(bas).modelMo = tempModelMoRes;
      
    calibStruct(bas).observed = OBS';
    calibStruct(bas).xlsResultsRow = num2str( (bas-1)*2+2);
    
    % ESTABLISH THE OBJECTIVE FUNCTION
    
    %To find month difference model v. observed
    diff  = calibStruct(bas).observed - calibStruct(bas).modelMo;
    diff2 = diff.^2;
    
    tssum   = sum(diff2);
    obj =  tssum;
    
    %display(obj);
    XX = sum(calibStruct(bas).modelMo);  %was  sum(calibStruct(bas).model)
    YY = sum(calibStruct(bas).observed);
    error = (YY-XX)/YY;
    
    GAUGE= calibStruct(bas).observed;
    MOD = ones(numMonthsRun,2); %BBB mod -- "numMonths" was previously 12
    MOD(:,2) = calibStruct(bas).modelMo;
    
    [b,bint,r,rint,stats] = regress(GAUGE,MOD);
    
    slope = b(2);
    r2 = stats(1);
    
    %Nash tells you how you did from a variance perspective
    nashdenom = var(GAUGE,1);
    nashnum = sum(diff2)/length(GAUGE);
    
    nash = (1 - (nashnum/nashdenom));
    
    %display(slope),display(r2),display('error(ob-wb)');display(error);
    
    calibStruct(bas).obj = obj;
    calibStruct(bas).slope = slope;
    calibStruct(bas).r2 = r2;
    calibStruct(bas).nash = nash;
    calibStruct(bas).error = error;
    
end

%msgbox('CALIBRATION COMPLETE.','INFO BOX');
