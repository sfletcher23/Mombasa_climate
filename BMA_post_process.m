%% BMA post processing


% NU and MU matrices are numSamples (10000) x numDecades (10) x numScenarios (9)

%% Trace plots

decades = {'1980-2000', '1990-2010', '2000-2020', '2010-2030', '2020-2040', ...
    '2030-2050', '2040-2060', '2050-2070', '2060-2080', '2070-2090'};

% Future precip
if true
for j = 1:9 % make a new figure for each scenario
    figure;
    for i = 1:10 % make a new subplot for each decade
        subplot(2,5,i)
        plot(1:1000,NU_futPrec(:,i,j))
        title(decades{i})
        ylim([0, 300])
        if i == 1 || i == 6
            ylabel('P (mm)')
        end
        if i >= 6
            xlabel('iteration')
        end
    end
    suptitle(strcat('Future P Trace Plot: Scenario:', ' ', num2str(j)))
end
end

% Future temp
if true
for j = 1:9 % make a new figure for each scenario
    figure;
    for i = 1:10 % make a new subplot for each decade
        subplot(2,5,i)
        plot(1:1000,NU_futTemp(:,i,j))
        title(decades{i})
        ylim([26, 31])
        if i == 1 || i == 6
            ylabel('T (degrees C)')
        end
        if i >= 6
            xlabel('iteration')
        end
    end
    suptitle(strcat('Future T Trace Plot: Scenario:', ' ', num2str(j)))
end
end

% Hist precip
if false
for j = 1:9 % make a new figure for each scenario
    figure;
    for i = 1:10 % make a new subplot for each decade
        subplot(2,5,i)
        plot(1:1000,MU_HistPrec(:,i,j))
        title(decades{i})
        ylim([0, 300])
        if i == 1 || i == 6
            ylabel('P (mm)')
        end
        if i >= 6
            xlabel('iteration')
        end
    end
    suptitle(strcat('Hist P Trace Plot: Scenario:', ' ', num2str(j)))
end
end

% Hist temp
if false
for j = 1:9 % make a new figure for each scenario
    figure;
    for i = 1:10 % make a new subplot for each decade
        subplot(2,5,i)
        plot(1:1000,MU_HistTemp(:,i,j))
        title(decades{i})
        ylim([20, 31])
        if i == 1 || i == 6
            ylabel('T (degrees C)')
        end
        if i >= 6
            xlabel('iteration')
        end
    end
    suptitle(strcat('Hist T Trace Plot: Scenario:', ' ', num2str(j)))
end
end

%% Trace plot zoom
% Comparing last decade with other decade at higher resolution to assess
% cycle problem

% Future Precip
sampleIndex = 900:1000; % Using a subset of the samples for better visualization
j = 1; % Scenario 
i = 10; % Decade 1
k = 9; % Deacde 2
figure;
subplot(2,1,1)
plot(sampleIndex, NU_futPrec(sampleIndex,i,j),'-o')
title(decades{i})
subplot(2,1,2)
plot(sampleIndex, NU_futPrec(sampleIndex,k,j),'-o')
title(decades{k})
suptitle(strcat('Zoomed trace plot fut P scenario:', num2str(j)))
% I played around with looking at different decades and scenarios - can
% really see the 21-step cycle. Although it doesn't draw the same exaclty values each time -
% suggests it is cycling through the models in a purely deterministic
% fashion but sampling with variation from each. Can we look at the marginal samples
% from another parameter, maybe lambda?, to confirm this? 

% Future Temp
sampleIndex = 900:1000; % Using a subset of the samples for better visualization
j = 1; % Scenario 5
figure;
subplot(2,1,1)
plot(sampleIndex, NU_futTemp(sampleIndex,i,j),'-o')
title(decades{i})
subplot(2,1,2)
plot(sampleIndex, NU_futTemp(sampleIndex,k,j),'-o')
title(decades{k})
suptitle(strcat('Zoomed trace plot fut T scenario 5:', num2str(j)))
% Can clearly see the difference between cycle and non-cycle


%% Autocorrelation

step = 1; % Increase the step to discard more of the samples
sampleIndex = 1:step:1000; 

% Future precip
for j = 1:9 % make a new figure for each scenario
    figure;
    for i = 1:10 % make a new subplot for each decade
        subplot(2,5,i)
        autocorr(NU_futPrec(sampleIndex,i,j),44)
        title(decades{i})
    end
    suptitle(strcat('Future P Autocorrelation: Scenario:', ' ', num2str(j)))
end

% Future temp
for j = 1:9 % make a new figure for each scenario
    figure;
    for i = 1:10 % make a new subplot for each decade
        subplot(2,5,i)
        autocorr(NU_futTemp(sampleIndex,i,j), 44)
        title(decades{i})
    end
    suptitle(strcat('Future T Autocorrelation: Scenario:', ' ', num2str(j)))
end




    


% 