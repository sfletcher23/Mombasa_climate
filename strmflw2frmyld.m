function  [safe_yield]  = strmflw2frmyld(inflow, storage)

maxRelease = 1000; % MCM/y , equal to 250,000 cm/d
for r = 40:maxRelease
    maxK = sequent_peak2(inflow, r);
    if maxK >= storage
        safe_yield = r;
        return
    end
end

error('didnt find safe yield')
    