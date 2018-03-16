%% Combine runoff

filenames = {'runoff_by_state_66414_15_Mar_2018_19_18_19', 'runoff_by_state_66413_15_Mar_2018_19_18_21'};

runoff_post = cell(M_T_abs, M_P_abs, N);
T_ts_post = cell(M_T_abs, N);
P_ts_post = cell(M_P_abs, N);

for k = 1:length(filenames)
    
    load(filenames{k},'runoff')
    
    for t = 1:N
        
        index_s_p_thisPeriod = index_s_p_time{t}; 
        for index_s_p = index_s_p_thisPeriod
            
            index_s_t_thisPeriod = index_s_t_time{t}; 
            for i= 1:length(index_s_t_thisPeriod)
                
                runoff_post{index_s_t_thisPeriod(i),index_s_p,t} = ...
                    [runoff_post{index_s_t_thisPeriod(i),index_s_p,t}; runoff{index_s_t_thisPeriod(i),index_s_p,t}];
                
            end
        end
    end
    
end
runoff = runoff_post;


for k = 1:length(filenames)
    
    load(filenames{k},'T_ts')
    
    for t = 1:N
            
        index_s_t_thisPeriod = index_s_t_time{t}; 
        for i= 1:length(index_s_t_thisPeriod)


            T_ts_post{index_s_t_thisPeriod(i),t} = ...
                [T_ts_post{index_s_t_thisPeriod(i),t}; T_ts{index_s_t_thisPeriod(i),t}];

        end

    end
    
end
T_ts = T_ts_post;



for k = 1:length(filenames)
    
    load(filenames{k}, 'P_ts')
    
    for t = 1:N
        
        index_s_p_thisPeriod = index_s_p_time{t}; 
        for index_s_p = index_s_p_thisPeriod

            P_ts_post{index_s_p,t} = ...
                [P_ts_post{index_s_p,t}; P_ts{index_s_p,t}];

        end
    end
    
end
P_ts = P_ts_post;



