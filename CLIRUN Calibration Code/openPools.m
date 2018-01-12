% openPools
isOpen = matlabpool('size') > 0;
while isOpen == 0
  matlabpool open 4
  isOpen = matlabpool('size') > 0;
end