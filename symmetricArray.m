function [ array, edgeStep ] = symmetricArray( stepSize, arraySize, column, offset )
% Creates an array centred symmetrically around zero, with arraySize number
% of elements and a given step size.
%
% The optional `column` argument selected whether the output array is a
% column (true) or row vector (false).
%
% If the `arraySize` argument is given as an N-element vector where N>1,
% then a cell array will be generated where each cell will be the array for
% the corresponding size argument. The `stepSize` and `offset` parameters
% can be 1 or N-element, where if 1-element same value will be used for all
% directions. To convert to 2D arrays, you can pass the cell arrays into
% meshgrid or similar.
%

% If not provided, assume an offset of zero
if (nargin < 4)
    offset = 0;
end

% Ensure column is a nice logical value
if (nargin > 2) && column
    column = true;
else
    column = false;
end

nDims = numel(arraySize);
if ~isscalar(stepSize) && (numel(stepSize) ~= nDims)
    error('Step size must be scalar or have one value per element in the `arraySize` argument');
end
if ~isscalar(offset) && (numel(offset) ~= nDims)
    error('Offset must be scalar or have one value per element in the `arraySize` argument');
end

% Determine the edge of the array, adding offset to both sides
baseEdgeStep = ((arraySize(:)-1)/2) .* stepSize(:);

% Calculate each array
array = cell(1,nDims);
edgeStep = zeros(1,nDims);
for dim = 1:nDims
    dimArray = makeArray(baseEdgeStep(dim), stepSize(min(end,dim)), arraySize(dim), offset(min(end,dim)));

    % Check if we need to transpose the array
    if column
        dimArray = dimArray.';
    end
    
    % Save edge step (after offset applied)
    edgeStep(dim) = dimArray(end);
    
    % Save to cell
    array{dim} = dimArray;
end

% For single dimension, maintain old behaviour or returning an array.
if isscalar(array)
    array = array{1};
end

end

function array = makeArray(edgeStep, stepSize, arraySize, offset)
% Generate an array for one direction.

% Calculate half the array
halfArray = -edgeStep : stepSize : 0;

% Then flip the array, accounting for odd/even length
if mod(arraySize,2)
    halfArray(end) = [];
    array = [halfArray 0 -fliplr(halfArray)];
else
    array = [halfArray -fliplr(halfArray)];
end

% Shift the array by offset number of steps
offset = offset * stepSize;
array = array + offset;

end

