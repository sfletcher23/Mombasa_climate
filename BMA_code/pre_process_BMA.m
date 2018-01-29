load('Input/Mombasa_TandP.mat'); 
rewrite_XYlambda = true;
create_scens = true;
reading_Routput = 0;
create_timeseries = 0;

%% Creating initial X, Y and lambda values.  These don't change 
if rewrite_XYlambda
    for year = 70:100
        YT0(year-69) = mean(T0(12*(year-1)+1:12*year));
        YP0(year-69) = log(mean(P0(12*(year-1)+1:12*year)));
    end

    X0 = [mean(YT0), mean(YP0)]';
    SD = [std(YT0), std(YP0)]';
    lambda0 = SD.^(-2);

    for year = 1:200
        YTij(year,:) = mean(Tij(12*(year-1)+1:12*year,:),1);
        YPij(year,:) = log(mean(Pij(12*(year-1)+1:12*year,:),1));
    end

    % This is when we use 1970-2000 to predict 2070-2100
    %X = [mean(YTij(70:100,:),1)', mean(YPij(70:100,:),1)']';
    %Y = [mean(YTij(FYi:FYi+30,:),1)', mean(YPij(FYi:FYi+30,:),1)']';

    % Saving variables for R code: 
    csvwrite('Input/lambda0.csv',lambda0)
    csvwrite('Input/Tij.csv',Tij)
    csvwrite('Input/Pij.csv',Pij)
    
    for decade = 1:11
        X = [mean(YTij(10*(decade-1)+80:10*(decade-1)+100,:),1)', mean(YPij(10*(decade-1)+80:10*(decade-1)+100,:),1)']';
        str1 = sprintf('Input/X_%2.0f.csv',1990+10*(decade-1));
        csvwrite(str1,X)
        %csvwrite(str2,Y)
    end
end

%% Creating scenarios
if create_scens
    deltaT_X0(1:10,1) = zeros(10,1);
    deltaT_X0(1:10,2) = 0:0.25:0.25*9;
    deltaT_X0(1:10,3) = 0:0.5:4.5;

    deltaP_X0(1:10,1) = 0:-0.08:-0.72;
    deltaP_X0(1:10,2) = zeros(10,1);
    deltaP_X0(1:10,3) = 0:0.08:0.72;

    % Univariate case
    X0_Topts = X0(1) + deltaT_X0;
    X0_Popts = log(exp(X0(2))+exp(X0(2))*deltaP_X0);

    csvwrite('Input/X0TU.csv',X0_Topts)
    csvwrite('Input/X0PU.csv',X0_Popts) 
end