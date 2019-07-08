clear all
close all

load('/Users/sarahfletcher/Dropbox (MIT)/Mombasa_Climate/Data/Mombasa_TandP.mat')

% Absolute change
figure;
subplot(2,2,1)
plot([1:2400]/12 + 1900,movmean(Pij_45,12*20))
title('RCP 4.5 Precp')
ylabel('mm/month')
ylim([40, 200])
subplot(2,2,2)
plot([1:2400]/12 + 1900,movmean(Pij_85,12*20))
title('RCP 8.5 Precp')
ylabel('mm/month')
ylim([40, 200])
subplot(2,2,3)
plot([1:2400]/12 + 1900,movmean(Tij_45,12*20))
title('RCP 4.5 Temp')
ylabel('^oC')
ylim([22, 30])
subplot(2,2,4)
plot([1:2400]/12 + 1900,movmean(Tij_85,12*20))
title('RCP 8.5 Temp')
ylabel('^oC')
ylim([22, 30])
suptitle('Absolute Comparison')

% Perc delta P and delta T
P0_45 = mean(Pij_45(1:20,:));
P0_85 = mean(Pij_85(1:20,:));
T0_45 = mean(Tij_45(1:20,:));
T0_85 = mean(Tij_85(1:20,:));
figure;
subplot(2,2,1)
plot([1:2400]/12 + 1900, (movmean(Pij_45,12*20) - P0_45) ./ P0_45)
title('RCP 4.5 Precp')
ylabel('Percent change from P0')
ylim([-.4, 1])
subplot(2,2,2)
plot([1:2400]/12 + 1900, (movmean(Pij_85,12*20) - P0_85) ./ P0_85)
title('RCP 8.5 Precp')
ylabel('Percent change from P0')
ylim([-.4, 1])
subplot(2,2,3)
plot([1:2400]/12 + 1900, movmean(Tij_45,12*20) - T0_45)
title('RCP 4.5 Temp')
ylabel('^oC')
ylim([-.5, 6])
subplot(2,2,4)
plot([1:2400]/12 + 1900, movmean(Tij_85,12*20) - T0_85)
title('RCP 8.5 Temp')
ylabel('^oC')
ylim([-.5, 6])
suptitle('Delta Comparison')