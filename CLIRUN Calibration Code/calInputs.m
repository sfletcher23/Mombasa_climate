%% Calibration inputs

% load temp, precip, trange, runoff
load([home.data '\Meteorological Data\CRU_Historical_Data\temp_ts2_1']);
load([home.data '\Meteorological Data\CRU_Historical_Data\trange_ts2_1']);
load([home.data '\Meteorological Data\CRU_Historical_Data\precip_ts2_1']);

load([home.data '\Runoff Data\GRDC\GRDC_CRU']);

temp2 = double(temp(crus, calMos)/10);
trange2 = double(trange(crus, calMos)/10);
precip2 = double(precip(crus, calMos)/10);

runoff = GRDC_CRU(crus, :);
runoff(runoff < 0) = NaN;
runoff(runoff > 0 & runoff < 1e-3) = NaN;
vec = mean(runoff, 2)./mean(precip2, 2);  %remove anomalous specific runoffs
runoff(vec < 0.1, :) = NaN;
runoff(vec == 1, :) = NaN;

tmeanBasins = extract_var_CRU(temp2, weights);
precipBasins = extract_var_CRU(precip2, weights);
trangeBasins = extract_var_CRU(trange2, weights);
runoffBasins = extract_var_CRU(runoff, weights);

% fill in missing basin 13 with neighbors
runoffBasins(10, :) = runoffBasins(9, :); 
runoffBasins(11, :) = mean([runoffBasins(9, :); runoffBasins(14, :)], 1);
runoffBasins(12, :) = runoffBasins(11, :);
runoffBasins(13, :) = mean([runoffBasins(14, :); runoffBasins(15, :)], 1);

clear temp temp2 precip precip2 trange trange2 GRDC_CRU runoff GRDCdata

%calculate PET
WATYEAR = 1;
petBasins = double(ModHargreaves3RSmat(lats,WATYEAR,tmeanBasins,...
  trangeBasins,precipBasins));

save([home.inputs '\MAT files\Madagascar34Basins_CalInputs'],...
  'tmeanBasins', 'precipBasins', 'trangeBasins', 'runoffBasins',...
  'petBasins', 'lats', 'weights', 'crus', 'startCalibYear', 'endCalibYear');

