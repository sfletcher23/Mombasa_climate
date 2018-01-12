%qcGraphics

bar(mean(runoffBasins, 2)./mean(precipBasins, 2))

% Maps comparing GRDC to precip
load([home.data '\Meteorological Data\CRU_Historical_Data\precip_ts2_1']);
mapCRU(mean(precip(:, 60*12+1:70*12), 2)/10, '-jet', [0 250]);

load([home.data '\Runoff Data\GRDC\GRDC_CRU']);
GRDC_CRU(GRDC_CRU < 0) = NaN;
mapCRU((nanmean(GRDC_CRU, 2) == 0)+1, '-jet');

%graphs comparing observed and modeled runoff

plot([calibStruct(i,1).observed calibStruct(i,1).model]);
if i == 2
  legend({'Obs','Mod'},'orientation','horiz','location','n');
end
xlim([1 84]);
ylabel(['Basin R',num2str(i)]);


%bars of model performance
temp.data = zeros(num.basins, 4);
temp.title = {'R2','MAR','Nash','Error'};
temp.data(i,1) = calibStruct(i,1).r2;
temp.data(i,2) = 12*mean(calibStruct(i,1).observed);
temp.data(i,3) = calibStruct(i,1).nash;
temp.data(i,4) = calibStruct(i,1).error;

for i = 1:4
    subplot(2,2,i)
    bar(temp.data(:, i));
    title(temp.title{i});
    xlim([0 num.basins+1]);
end
pet_day = PET(1:84)./repmat(days, 1, 7);
