%% Post processing BMA data

% Which files to open - need to set this every run
dateOpen = '2018-01-28';
% Need to input jobid of T and P runs if running on cluster
if exist(getenv('SLURM_JOB_ID'))
    jobIdT = 12345;
    jobIdP = 12345;
else
    jobIdT = '';
    jobIdP = '';
end

reading_Routput = true;
create_timeseries = true;
make_plots = false;

%% Reading R Output and creating time series 
if reading_Routput

    for scen_ii = 1%:3
%         filepath = scen_dir_names{scen_ii};
%         if ~exist(filepath,'dir') %if doesn't exist
%             mkdir(filepath)
%         end
%         close all
        time_ii = 0;
        for year = 2000%:10:2090
            time_ii = time_ii+1;
            tmpstr = strcat('Output/', sprintf('muUT_%d_scen%d',year,scen_ii),'_','job', '_',jobIdT,'_', dateOpen,'_.csv');
            tmp = csvread(tmpstr);
            MUT(:,time_ii,scen_ii) = tmp;

            tmpstr = strcat('Output/',sprintf('mnUT_%d_scen%d',year,scen_ii),'_','job', '_',jobIdT,'_', dateOpen,'_.csv');
            tmp = csvread(tmpstr); 
            NUT(:,time_ii,scen_ii) = tmp; 

            tmpstr = strcat('Output/',sprintf('muUP_%d_scen%d',year,scen_ii),'_','job', '_',jobIdP,'_', dateOpen,'_.csv');
            tmp = csvread(tmpstr);
            MUP(:,time_ii,scen_ii) = exp(tmp);

            tmpstr = strcat('Output/',sprintf('nuUP_%d_scen%d',year,scen_ii),'_','job', '_',jobIdP,'_', dateOpen,'_.csv');
            tmp = csvread(tmpstr); 
            NUP(:,time_ii,scen_ii) = exp(tmp); 
        end
    end
    
    saveName = strcat('BMA_results', datestr(datetime)); 
    save(saveName,'MUP','MUT','NUT','NUP', 'jobIdT', 'jobIdP')
    
end

%% Make plots

% Make sure plotting off if running on svante to avoid error
if exist(getenv('SLURM_JOB_ID'))
    make_plots = false;
end

if make_plots

    scen_titleT{1} = 'scenario 1-T, \Delta T = 0';
    scen_titleP{1} = 'scenario 1-P, \DeltaP = -72%';
    scen_titleT{2} = 'scenario 2-T, \Delta T = 2.25';
    scen_titleP{2} = 'scenario 2-P, \DeltaP = 0%';
    scen_titleT{3} = 'scenario 3-T, \Delta T = 4.5';
    scen_titleP{3} = 'scenario 3-P, \DeltaP = +72%';

    for scen_ii = 1:3
        coloropts = [1,6,11,16,21,36,41,46,51]; % with time_jj up to 9
        temp_opts1 = [24,24,24];
        fig = figure
            str = sprintf('End2080_scen%d.tif',scen_ii);
            baseFileName = sprintf(str);
            fullFileName = fullfile(filepath,baseFileName); 
            set(gcf, 'PaperUnits', 'centimeters');
            set(gcf, 'Position', [0 0 700 700])
        subplot(2,1,1)
        for time_jj = 1:9
            histf(NUT(:,time_jj,scen_ii),temp_opts1(scen_ii):.03:temp_opts1(scen_ii)+6,'facecolor',ColorMeg(coloropts(time_jj),:),'facealpha',.5,'edgecolor',[0.3,0.3,0.3])
            hold on
        end
        alpha(0.9); box on; axis tight; xlabel('annual temperature(^oC)'); ylabel('frequency (in 1000)')
        legend('1980-2000','1990-2010','2000-2020','2010-2030','2020-2040','2030-2050','2040-2060','2050-2070','2060-2080','location','northwest')
        legend boxoff; title(scen_titleT{scen_ii});

        subplot(2,1,2)
        for time_jj = 1:9
            histf(NUP(:,time_jj,scen_ii),50:0.5:90,'facecolor',ColorMeg(coloropts(time_jj),:),'facealpha',.5,'edgecolor',[0.3,0.3,0.3])
            hold on
        end
        alpha(0.9); box on; axis tight; xlabel('annual precip (mm/month)'); ylabel('frequency (in 1000)')
        set(fig, 'PaperPositionMode', 'auto');
        title(scen_titleP{scen_ii});

    end
    
end


%% Timeseries - need to read in initial data
% creating standardized anomalies relative to decadal average 
    % (From Sarah: I'm sure there is a way to vectorize this that would be
    % more efficient, but this is definitely a small piece of the overall
    % cpu time so prob not worth worrying 
if create_timeseries
    for decade = 1:20
        for model = 1:21    
            Tmean = mean(Tij(120*(decade-1)+1:120*decade,model)); %decadal mean
            Pmean = mean(Pij(120*(decade-1)+1:120*decade,model));
            % anomalies relative to decadal mean
            T_anoms(120*(decade-1)+1:120*decade,model) = Tij(120*(decade-1)+1:120*decade,model)-Tmean;
            P_anoms(120*(decade-1)+1:120*decade,model) = (Pij(120*(decade-1)+1:120*decade,model)-Pmean)/Pmean;    
        end
    end


    for decade = 1:10
        for scen_tt = 1:3
            for scen_pp = 1:3
                %for each combination of tt and pp scenario, we choose 100 of
                %the 1000 draws from NUT and NUP
                r = randi([1 1000],1,100);
                for drawi = 1:100
                    for year = 1:10
                        m = randi([1 21]); % choosing model number
                        y = randi([1 10]);
                        Tmean = NUT(r(drawi),decade,scen_tt);
                        oind = 12*(year-1)+1;
                        iind = 1200+120*(decade-1)+12*(y-1)+1; % Here, we skip the first 100 years of anomalies
                        T_timeseries(oind:oind + 11,drawi,decade,3*(scen_pp-1)+scen_tt) = T_anoms(iind:iind+11,m)+Tmean;
                        Pmean = NUP(r(drawi),decade,scen_pp);
                        P_timeseries(oind:oind + 11,drawi,decade,3*(scen_pp-1)+scen_tt) = Pmean*P_anoms(iind:iind+11,m)+Pmean;
                    end
                end
            end
        end
    end
    save(saveName,'T_timeseries','P_timeseries')
end