function  [maxK, K]  = sequent_peak2(inflow, release)

% Inflow is time series and Release is constant with same volumetric units

K = zeros(1,length(inflow)+1);
K0 = 0;
for t = 1:length(inflow)
    if t == 1
        K(t) = max(0,K0 + release(t) - inflow(t));
    else
    K(t) = max(0,K(t-1) + release(t) - inflow(t));
    end
end
maxK = max(K);