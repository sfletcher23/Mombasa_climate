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
    indexFlex = find(costmodel.storage == flex_storage);
    ecost = added_storage * costmodel.unit_cost(indexFlex)*1.5;
    dcost = dcost;
end


end