function [elementPositions, cncPositions] = removeNaN(elementPositions, cncPositions)
% Removes NaN values from element/CNC positions arrays.

nCncDims = size(cncPositions,3);

cncPositions(repmat(isnan(elementPositions),1,1,nCncDims)) = [];
cncPositions = reshape(cncPositions,[],nCncDims);

elementPositions(isnan(elementPositions)) = [];
elementPositions = reshape(elementPositions, [], 1);

end
