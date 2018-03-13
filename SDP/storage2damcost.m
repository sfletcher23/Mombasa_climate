function [dcost, ecost] = storage2damcost(storage, flex_storage)

load('dam_cost_model')
ecost = [];

if ~ismember(storage, costmodel.storage)
   error('invaild storage volume')
end

index = find(costmodel.storage == storage);
dcost = costmodel.dam_cost(index);

if flex_storage
    added_storage = flex_storage - storage;
    ecost = added_storage * costmodel.exp_unit_cost(index);
    dcost = dcost * 1.1;
end


end