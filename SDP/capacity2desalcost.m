function [dcost, ecost, opex] = capacity2desalcost(capacity, flex_capacity)

ecost = [];

% capex
unitcost = -189.1 * log(mcmpy2cmpd(capacity)) + 3385;
dcost = unitcost * mcmpy2cmpd(capacity);

% opex
opex = -0.0001 * mcmpy2cmpd(capacity) + 186.76;

if flex_capacity
    unitcost = -189.1 * log(mcmpy2cmpd(flex_capacity)) + 3385;
    dcost = unitcost * mcmpy2cmpd(flex_capacity);
    totalcost = dcost * 1.05;   % 5% more expensive; compared to 8% for Melbourne going from 75 to 150
    ecost = totalcost * .2; 
    dcost = totalcost * .8;
    opex = -0.0001 * mcmpy2cmpd(flex_capacity) + 186.76;
end

end