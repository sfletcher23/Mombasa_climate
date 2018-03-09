%% Combine runoff

filenames = {'runoff_by_state_63724_01_Mar_2018_12_00_34', 'runoff_by_state_63681_28_Feb_2018_09_06_37', ...
    'runoff_by_state_63728_01_Mar_2018_12_01_33', 'runoff_by_state_63725_01_Mar_2018_14_03_01',...
    'runoff_by_state_63726_01_Mar_2018_14_00_25', 'runoff_by_state_63727_01_Mar_2018_13_56_16',...
    'runoff_by_state_63791_01_Mar_2018_22_11_32'};

runoff_post = cell(M_T_abs, M_P_abs, N);
for k = 1:length(filenames)
    
    load(filenames{k},'runoff')
    
    for t = 1:N
        
        index_s_p_thisPeriod = index_s_p_time{t}; 
        for index_s_p = index_s_p_thisPeriod
            
            index_s_t_thisPeriod = index_s_t_time{t}; 
            for i= 1:length(index_s_t_thisPeriod)
                
                runoff_post{index_s_t_thisPeriod(i),index_s_p,t} = ...
                    [runoff_post{index_s_t_thisPeriod(i),index_s_p,t}; runoff{i,index_s_p,t}];
           
            end
        end
    end
    
end
runoff = runoff_post;


