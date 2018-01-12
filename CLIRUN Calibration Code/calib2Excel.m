% calib2Excel

tic
basinNames = basins';
num.basins = length(basins);

BASIN = fullfile(home.outputs,'Runoff', 'calibResults.xlsx');
header1 = {'basin' 't_low' 't_high' 'ku' 'kp' 'kl' 'sat' 'interc' 'over' 'd melt'};
xlswrite(BASIN, header1,'calibration', 'a1:j1');
header2 = {'obj' 'slope' 'r2' 'Nash-S' 'error(ob-wb)'};
xlswrite(BASIN, header2,'calibration','k1:o1');
header3 = {'basin' 'jan' 'feb' 'mar' 'apr' 'may' 'jun' 'jul' 'aug'...
    'sep' 'oct' 'nov' 'dec'};
xlswrite(BASIN, header3,'results','a1');
header4 = {};
header5 = {};
for bas = 1:length(basins)
    iStr = num2str(basins(bas));
    header4 = vertcat( header4 ,{['mod-',iStr] ;['obs-',iStr] } );
    header5 = vertcat( header5 ,{['modRun-', iStr] ;['obsPrecip-',iStr] } );
end
xlswrite(BASIN, header4,'results','a2');
xlswrite(BASIN, header5,'resultsDay','a2');
xlswrite(BASIN, basinNames, 'calibration', 'a2');       



try 
  temp.obs = zeros(num.basins, length(calibStruct(1).modelMo));
  ind1 = true;
  lenM = length(calibStruct(1).modelMo);
catch
  temp.obs = zeros(num.basins, length(calibStruct(1).model));
  ind1 = false;
  lenM = length(calibStruct(1).model);
end
temp.mod = temp.obs; 
if ind1
  temp.modDay = zeros(num.basins, length(calibStruct(1).model));
end
temp.genes = zeros(num.basins, 9);
temp.perform = zeros(num.basins, 5);
for bas=[1:length(basins)]
  if ind1
    temp.mod(bas, :) = calibStruct(bas).modelMo';
  else temp.mod(bas, :) = calibStruct(bas).model';
  end
  temp.obs(bas, :) = calibStruct(bas).observed';
  if ind1
    temp.modDay(bas, :) = calibStruct(bas).model';
  end
  temp.genes(bas, 1:9) = calibStruct(bas).x;
  temp.perform(bas, 1)= calibStruct(bas).obj;
  temp.perform(bas, 2)= calibStruct(bas).slope;
  temp.perform(bas, 3)= calibStruct(bas).r2;
  temp.perform(bas, 4)= calibStruct(bas).nash;
  temp.perform(bas, 5)= calibStruct(bas).error;
end

temp.results = reshape(permute(cat(3, temp.mod, temp.obs), [3 1 2]),...
  [], lenM);
if ind1
  temp.resultsDay = reshape(permute(cat(3, temp.modDay, reshape(precip, [], length(basins))'),...
    [3 1 2]), [], length(calibStruct(1).model));
end
xlswrite(BASIN, [temp.genes temp.perform],'calibration', 'b2');
xlswrite(BASIN, temp.results,'results', 'b2');
if ind1
  xlswrite(BASIN, temp.resultsDay,'resultsDay', 'b2');
end


disp(['Writing calib results to excel took: ',num2str(toc) ' seconds']);
