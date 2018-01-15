%% Hydrological data processing

% Streamflow is a numYears(14) x numMonths(12) matrix with monthly streamflow data
% in m^3/s. The data ranges from June 1976 to May 1990. The first column is
% June; the last column is May. 


% P0 contains monthly historical precipitation data in mm/month at Mombasa from CRU
% from 1901 - 2012

% T0 contains monthly historical temperature data in degrees C at Mombasa
% from CRU from 1901 - 2012

%% Load data
load('historical data')

%% Select matching date range from CRU data

P0_month = reshape(P0, 12, [])'; 
index76 = 1976 - 1901 + 1; 
index90 = 1990 - 1901 + 1; 
P0_76_to_90 = P0_month(index76:index90,:);

datestart = (1976-1901)*12 + 6; 
dateend = (1990-1901)*12 + 5;

P0_date_lin = P0(datestart:dateend);
P0_date = reshape(P0_date_lin, 12, [])';

T0_date_lin = T0(datestart:dateend);
T0_date =  reshape(T0_date_lin, 12, [])';

Tmin_date_lin = Tmin(datestart:dateend);
Tmin_date =  reshape(Tmin_date_lin, 12, [])';

Tmax_date_lin = Tmax(datestart:dateend);
Tmax_date =  reshape(Tmax_date_lin, 12, [])';

PET_date_lin = PET(datestart:dateend);
PET_date = reshape(PET_date_lin, 12, [])';

%% Calcuate PET using Modified Hargreaves

lat = -3.96;   % Approx lat at proposed dam site
watyear = 6; % Data starts in June
temp = T0_date;
temprange = Tmax_date - Tmin_date;
prec = P0_date;
pet_hg = ModHargreaves4(ones(14,1)*lat,watyear,temp,temprange,prec);

%% Save pre processed data for CLIRUN
area = 2250 * 1E6; %m2
runoff = reshape(streamflow',1,[]) * 60 * 60 * 24 * 365  / area; % converting from m3/s to m3/y and dividing by area
runoff = runoff * 1E3; % Converting from m/y to mm/y;
save('Mwache_hydro_data', 'PET_date_lin', 'P0_date_lin', 'watyear', 'runoff')

%% Calculate storage-yield curve using sequent-peak algorithm

inflow = reshape(streamflow',1,[]);
%inflow = streamflow;
[numRuns,numTime] = size(inflow);
maxRelease = 9.448E6;
release = 0:maxRelease/25:maxRelease;
maxK = zeros(numRuns,length(release));
mar = zeros(numTime,1);
for i = 1:length(release)
    [maxK(:,i), mar]  = sequent_peak(inflow, release(i), 6);
end

if false
% Calculate inverse of sequent-peak
syms  x
sqpk = @(x) sequent_peak(inflow, x, 6);
storage = 200E6;
eqn = storage == sqpk(x);
solx = solve(eqn, x);
end
    
%% Plots

if true
% Historical P, T, streamflow
figure
subplot(3,2,1)
plot(1:12, streamflow)
xlim([1 12])
ax = gca;
ax.XTick = 1:12;
ax.XTickLabels = {'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M'};
title('Streamflow')
ylabel('cm/s')

subplot(3,2,2)
plot(1:12*14, reshape(streamflow, 1, []))
xlim([1 12*14])
ax = gca;
ax.XTick = 8:12:12*14;
ax.XTickLabels = num2cell(1977:1990);
title('Streamflow')
ylabel('cm/s')

subplot(3,2,3)
datestart = (1977-1901)*12+1; 
dateend = (1990-1901)*12 + 12;
plot(1:12, P0_date)
xlim([1 12])
ax = gca;
ax.XTick = 1:12;
ax.XTickLabels = {'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M'};
title('Precipitation')
ylabel('mm/month')

subplot(3,2,4)
plot(1:12*14, P0_date_lin)
xlim([1 12*14])
ax = gca;
ax.XTick = 8:12:12*14;
ax.XTickLabels = num2cell(1977:1990);
title('Precipitation')
ylabel('mm/month')

subplot(3,2,5)
plot(1:12,T0_date)
xlim([1 12])
ax = gca;
ax.XTick = 1:12;
ax.XTickLabels = {'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M'};
title('Temperature')
ylabel('degrees C')

subplot(3,2,6)
plot(1:12*14, T0_date_lin)
xlim([1 12*14])
ax = gca;
ax.XTick = 8:12:12*14;
ax.XTickLabels = num2cell(1977:1990);
title('Temperature')
ylabel('degrees C')

% Plot streamflow and P together monthly, each year in differnet panel
figure
for i = 1:14
    subplot(7, 2, i)
    yyaxis left
    plot(streamflow(i,:));
    ylim([0 20])
    hold on
    yyaxis right
    plot(P0_date(i,:));
    xlim([1 12])
    ylim([0 400])
    ax = gca;
    ax.XTick = 1:12;
    ax.XTickLabels = {'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M'};
end
legend('S','P')
suptitle('Comparing Streamflow and P')

% Plot Tmin and Tmax
figure
for i = 1:14
    subplot(7, 2, i)
    plot(Tmin_date(i,:));
    hold on
    plot(Tmax_date(i,:));
    plot(T0_date(1,:));
    ylim([10 max(max(Tmax_date))+5])
    xlim([1 12])
    ax = gca;
    ax.XTick = 1:12;
    ax.XTickLabels = {'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M'};
end
legend('Tmin','Tmax', 'Tavg')

% Compare PET from modHG and CRU data
figure;
for i = 1:14
    subplot(7, 2, i)
    plot(pet_hg(i,:));
    hold on
    plot(PET_date(i,:));
    xlim([1 12])
    ax = gca;
    ax.XTick = 1:12;
    ax.XTickLabels = {'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M'};
end
legend('PET HG', 'PET CRU')
suptitle('Comparing PET modeled estimates with CRU data')
er = pet_hg - PET_date;
sq_er = er .^ 2;
mean_sq_er = mean(mean(sq_er));
rmse = sqrt(mean_sq_er);
end

% Plot storage yield curve
figure
plot(maxK, release)
hold on
line([0 max(max(maxK))], [mar mar])
ylabel('Yield')
xlabel('Storage')

