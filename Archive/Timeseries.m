
% Timeseries for flow model

% creating standardized anomalies relative to decadal average
load('/Users/sarahfletcher/Documents/MATLAB/Mombasa_Climate/BMA_code/Input/Mombasa_TandP.mat')
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
    %save('TP_time_series_Jan27.mat','T_timeseries','P_timeseries')
end