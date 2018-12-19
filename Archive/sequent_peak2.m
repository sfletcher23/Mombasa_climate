function  [maxK, K]  = sequent_peak2(inflow, release)

% Inflow is time series and Release is constant with same volumetric units
if length(release) == 1
    release = release * ones(1,length(inflow));
end

K = zeros(1,length(inflow));
K0 = 0;
for t = 1:length(inflow)
    if t == 1
       % K(t) = feval(symengine,'max',0 - K0 + inflow(t),release); 
       K(t) = max(0,K0 + release(t) - inflow(t));
    else
        %K(t) = feval(symengine,'max',0, K(t-1) + release - inflow(t)); 
        K(t) = max(0,K(t-1) + release(t) - inflow(t));
    end
end
% sortK = sort(K);
% maxK = sortK(end);
maxK = max(K);