%% Hydrological data processing

% Streamflow is a numYears(14) x numMonths(12) matrix with monthly streamflow data
% in m^3/s. The data ranges from June 1976 to May 1990. The first column is
% June; the last column is May. 


% P0 contains monthly historical precipitation data in mm/month at Mombasa from CRU
% from 1901 - 2012

% T0 contains monthly historical temperature data in degrees C at Mombasa
% from CRU from 1901 - 2012

load('historical hydrological data')

%% Plot historical data

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
plot(1:12, reshape(P0(datestart:dateend),12,14)')
xlim([1 12])
ax = gca;
ax.XTick = 1:12;
ax.XTickLabels = {'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M'};
title('Precipitation')
ylabel('mm/month')

subplot(3,2,4)
plot(1:12*14, P0(datestart:dateend))
xlim([1 12*14])
ax = gca;
ax.XTick = 8:12:12*14;
ax.XTickLabels = num2cell(1977:1990);
title('Precipitation')
ylabel('mm/month')

subplot(3,2,5)
plot(1:12, reshape(T0(datestart:dateend),12,14)')
xlim([1 12])
ax = gca;
ax.XTick = 1:12;
ax.XTickLabels = {'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J', 'F', 'M', 'A', 'M'};
title('Temperature')
ylabel('degrees C')

subplot(3,2,6)
plot(1:12*14, T0(datestart:dateend))
xlim([1 12*14])
ax = gca;
ax.XTick = 8:12:12*14;
ax.XTickLabels = num2cell(1977:1990);
title('Temperature')
ylabel('degrees C')



