function  [yield, K, demand, unmet_dom, unmet_ag]  = runoff2yield(inflow, T, P, storage, runParam, climParam)

% Inflow is a monthly time series in MCM/y starting in January
% Storage is a scalar 

numYears = runParam.steplen;
[numRuns,~] = size(T);

% check numRuns correct
dmd_dom = cmpd2mcmpy(186000) * ones(numRuns,12*numYears);
dmd_dom = cmpd2mcmpy(150000) * ones(numRuns,12*numYears);
if runParam.desalOn
    dmd_dom = cmpd2mcmpy(400000);
end
dmd_ag = repmat([2.5 1.5 0.8 2.0 1.9 2.9 3.6 0.6 0.5 0.3 0.2 3.1], numRuns,numYears);
demand = dmd_dom + dmd_ag;
dead_storage = 20;
env_flow = 0;
eff_storage = storage - dead_storage;

[E]  = evaporation_sdp(storage, T, P, climParam, runParam);
    
K = zeros(numRuns,numYears*12);
release = zeros(numRuns,numYears*12);
K0 = eff_storage;

if ~runParam.desalOn
    for t = 1:numYears*12
        if t == 1
            Kprev = K0;
        else
            Kprev = K(t-1);
        end

        % If demand is less than effective inflow, release all demand and add storage up to limit
        indLess = demand(:,t) < inflow(:,t) - env_flow - E(:,t);
        release(indLess,t) = demand(indLess,t);
        K(indLess,t) = min(Kprev - release(indLess,t) + inflow(indLess,t) - env_flow - E(indLess,t), eff_storage);
        % If demand is greater than effective inflow, but less than available storage, release all demand    
        indMid = demand(:,t) < Kprev + inflow(:,t) - env_flow - E(:,t) & demand(:,t) > inflow(:,t) - env_flow - E(:,t);
        release(indMid,t) = demand(indMid,t);
        K(indMid,t) = Kprev - release(indMid,t) + inflow(indMid,t) - env_flow - E(indMid,t); 
        % If demand is greater than effective inflow and storage, release as much as available
        indGreat = ~indLess & ~indMid;
        release(indGreat,t) = Kprev + inflow(indGreat,t) - env_flow - E(indGreat,t);
        K(indGreat,t) = 0;

    end
else
    
desalsupply = zeros(numRuns,numYears*12);
desalthreshold = runParam.desalcapacity 
    for t = 1:numYears*12
        if t == 1
            Kprev = K0;
        else
            Kprev = K(t-1);
        end

        % If demand is less than effective inflow, release all demand and add storage up to limit
        indLess = demand(:,t) < inflow(:,t) - env_flow - E(:,t);
        release(indLess,t) = demand(indLess,t);
        K(indLess,t) = min(Kprev - release(indLess,t) + inflow(indLess,t) - env_flow - E(indLess,t), eff_storage);
        % If demand is greater than effective inflow, but less than available storage, release all demand    
        indMid = demand(:,t) < Kprev + inflow(:,t) - env_flow - E(:,t) & demand(:,t) > inflow(:,t) - env_flow - E(:,t);
        release(indMid,t) = demand(indMid,t);
        K(indMid,t) = Kprev - release(indMid,t) + inflow(indMid,t) - env_flow - E(indMid,t); 
        % If demand is greater than effective inflow and storage, release as much as available
        indGreat = ~indLess & ~indMid;
        release(indGreat,t) = Kprev + inflow(indGreat,t) - env_flow - E(indGreat,t);
        K(indGreat,t) = 0;
        
        

    end
end
    

% Ag demand is unmet first
unmet = max(demand - release, 0);
unmet_ag = min(unmet, dmd_ag);
unmet_dom = unmet - unmet_ag;

yield = release;

