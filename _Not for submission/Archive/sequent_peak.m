function  [maxK, MAR]  = sequent_peak(inflow, release, monthstart)

% Inputs:
% Inflows: monthly average streamflow in m3/sec  [numRuns x numTime] 
% Release: daily release [1 x 1] in m3/day
% monthstart: starting month eg June is 6

% Output:
% maxK  in m3
% MAR in m3/d


% number of time series input
[numRuns, numTime] = size(inflow);

% days in each month
days = [31 28 31 30 31 30 31 31 30 31 30 31]; % original
days = [days(monthstart:end) days(1:monthstart-1)]; % data starts in June
days = [repmat(days,1, floor(numTime/12))  days(1:mod(numTime,12))]; % for length of timeseries

% Inflow 
Q = inflow * 60 * 60 * 24 .* days; % m3/m

% MAR
MAR = mean(Q,2) / mean(days); % m3/d

% Release
R = release * days; % Release m3/m - from dam design spec in WB report

% Apply algorithm
K0 = 0;
K = zeros(numRuns,numTime);
for i = 1:numRuns
    for t = 1:numTime
        if t == 1
            K(i,t) = max(0, K0 + R(t) - Q(i,t));
        else
            K(i,t) = max(0, K(t-1) + R(t) - Q(i,t));
        end
    end
end

% symmax =  @(x,y,z)feval(symengine,'max',x,y,z);
% maxK = symmax(K,[],2);
maxK = max(K,[],2); % m3/m

