classdef GridPath < double
    % Utility for converting one or more arrays of points (one per dim)
    % into a single fully populated array of all points.

    enumeration
        Sparse      (0)  % Points simply used directly, each dim must be the same length
        MeshGrid    (1)  % Treated as regular grid for each dimension and formed into a mesh grid
        Diamond     (2)  % Treated as regular grid for each dimension and formed into a mesh grid, but keeps every other point in a diamond formation
        ZigZag      (3)  % Formed into a continuous zig-zag mesh grid, minimising changes when stepping each point
       %Hilbert     (4)  % Formed into a hilbert space filling curve
        SpiralGrid  (5)  % Spiral on grid in X-Y. For >2 dims, will switch to zig-zag for each extra dim.
        SpiralSlice (6)  % Sparsely forms a spiral in X-Y, X-Z then Y-Z planes which intersect in the centre.
    end


    methods (Static = true)
        % Generate X-Y-... points Helpers
        %  - For these inputs are vectors one for each dimension. These
        %    must be scalar (same for all dims) or at least nargout entries
        %  - Points for each dimension are returned as each output argument

        % Uniform grid from start to no more than stop in step increments
        function [varargout] = createUniformPointsRange(start, step, stop)
            nDims = nargout;
            varargout = cell(1,nDims);
            % Check inputs
            start = checkArgSize(start, nDims, 'start');
            step = checkArgSize(step, nDims, 'step');
            stop = checkArgSize(stop, nDims, 'stop');
            % Make outputs
            for dim = 1:nDims
                varargout{dim} = start(dim):step(dim):stop(dim);
            end
        end

        % Uniform grid from start in nSteps step-sized increments
        function [varargout] = createUniformPointsNSteps(start, step, nSteps)
            nDims = nargout;
            varargout = cell(1,nDims);
            % Check inputs
            start = checkArgSize(start, nDims, 'start');
            step = checkArgSize(step, nDims, 'step');
            nSteps = checkArgSize(nSteps, nDims, 'nSteps');
            % Make outputs
            for dim = 1:nDims
                varargout{dim} = start(dim) + step(dim).*(0:(nSteps(dim)-1));
            end
        end
        
        % Uniform grid out from centre in nSteps step-sized increments
        function [varargout] = createUniformPointsCentred(centre, step, nSteps)
            nDims = nargout;
            % Check inputs
            centre = checkArgSize(centre, nDims, 'centre');
            step = checkArgSize(step, nDims, 'step');
            nSteps = checkArgSize(nSteps, nDims, 'nSteps');
            % Make outputs
            varargout = symmetricArray(step, nSteps, false, centre ./ step);
        end

        % Generate Grid Helpers

        % Create the grid from points
        %  - varargin is each subsequent dimension if more than two are
        %    needed
        function [grid, nPoints, order] = createGrid(objOrGridNum, xPoints, yPoints, varargin)
            narginchk(3,inf);
            obj = GridPath.lookupObjOrNum(objOrGridNum);
            [nPoints, ~, grid, order] = obj.(['create' char(obj)])([], xPoints, yPoints, varargin{:});
        end

        % Check the size of the grid without actually making it
        function nPoints = gridCount(objOrGridNum, xPoints, yPoints, varargin)
            narginchk(3,inf);
            obj = GridPath.lookupObjOrNum(objOrGridNum);
            nPoints = obj.(['create' char(obj)])([], xPoints, yPoints, varargin{:});
        end
        
        % Select 2D slice best encapsulating grid index
        %  - This is intended for use with plotting, for example should
        %    something operating on a grid path plot say an X-Y slice (and if
        %    so which), or say an X-Z slice.
        %  - sliceDims is a 2-element array indicating which dimensions are
        %    used for the 2D slice
        %  - sliceColon is a cell array containing ':' for the selected
        %    slice dimensions, and scalar array indicies for the other
        %    dimensions to allow selecting the correct data subset.
        %  - gridIdx is the index of the point to be plotted. This can
        %    either be a scalar integer index, in which case grid and order
        %    will be calculated automatically. Alternatively it can be a 2
        %    element cell array of {index, order} to save time by not
        %    having to recalculate the grid order.
        function [sliceColon, sliceDims] = select2DSlice(objOrGridNum, gridIdx, xPoints, yPoints, varargin)
            narginchk(4,inf);
            obj = GridPath.lookupObjOrNum(objOrGridNum);
            % Find our slice
            [~,sliceColon] = obj.(['create' char(obj)])(gridIdx, xPoints, yPoints, varargin{:});
            % Grab the two slice dims from the colon array - there should
            % be exactly two char ':' in the colon cell.
            sliceDims = find(cellfun(@ischar, sliceColon));
        end
        
    end

    methods (Access = protected, Static = true)
        % Find enum object from enum, grid number, or char name.
        function obj = lookupObjOrNum(objOrGridNum)
            if isa(objOrGridNum,mfilename('class'))
                obj = objOrGridNum;
            else
                obj = feval(mfilename('class'), objOrGridNum);
            end
        end

        % The following are creators for each grid type. The trace of the
        % function must be:
        %
        %   [nPoints, slice, grid, order] = create<EnumName>(gridIdx, varargin)
        % 
        % Arguments:
        %
        %   - The gridIdx is used only for optionally generating the slice
        %   output.
        %
        %   - Subsequent varargins are a set of numeric points arrays, one
        %   for each dimension.
        %
        % Returns:
        %   - nPoints output indicates the number of points in the grid.
        %     When nargout == 1, try to optimise if possible to only giving
        %     the number of expected points and not generating the grid.
        %
        %   - grid is an [N x M] array where N is number of points and M is
        %     the number of dimensions set by input arguments
        %   
        %   - order is a map from current grid index to the corresponding
        %     point in a uniform grid of points. For example if accessing
        %     the gridIdx'th point of order gives the array index into the
        %     equivalent fully populated N-D array.
        %
        %   - slice is a cell array containing exactly two ':' strings, and
        %     nDims-2 scalar integer cells. This is used to select the
        %     appropriate data slice for plotting for a given grid index.
        %
        
        % Sparse creator
        function [nPoints, slice, grid, order] = createSparse(~, varargin)
            gridSize = meshGridSize(varargin{:});
            nPoints = max(gridSize);
            if nargout >= 3
                if ~all(gridSize == nPoints)
                    error('GridPath:InvalidPoints', 'For Sparse grids, all points arrays must have the same number of elements');
                end
                grid = flattenMerge(varargin);
                % Sparse is simply ordered exactly as given by the user
                order = columnColon(1,nPoints);
            end
            if nargout >= 2
                % Slice is basically everything as it has no real meaning
                % for sparse
                slice = defaultSlice(2);
            end
        end

        % MeshGrid creator
        function [nPoints, slice, grid, order] = createMeshGrid(gridIdx, x, y, varargin)
            gridSize = meshGridSize(x, y, varargin{:});
            nPoints = prod(gridSize);
            if nargout >= 3
                grids = cell(1, nargin-1);
                [grids{:}] = ndgrid(x, y, varargin{:});
                grid = flattenMerge(grids);
                % Mesh grids are simply in order.
                order = columnColon(1,nPoints);
            end
            if nargout >= 2
                % Mesh grids always select an xy slice.
                slice = xySlice(gridSize, gridIdx);
            end
        end

        % ZigZag creator
        function [nPoints, slice, grid, order] = createZigZag(gridIdx, x, varargin)
            if nargin == 2
                nPoints = numel(x);
                grid = x(:);
                order = columnColon(1,nPoints);
                slice = defaultSlice(2);
                return
            end
            gridSize = meshGridSize(x, varargin{:});
            nPoints = prod(gridSize);
            if nargout >= (2 + iscell(gridIdx))
                % Get the sub-grid for all remaining dims
                [~, ~, subGrid, subOrder] = GridPath.createZigZag(gridIdx, x, varargin{1:end-1});
                % Loop through all indicies in current dimension
                dimPoints = varargin{end};
                grid = [];
                order = [];
                for dimIdx = 1:numel(dimPoints)
                    % For each, add a new copy of the sub-grid for this new
                    % point. We flip the sub-grid for even points.
                    [newGrid, newOrder] = addNewDim(subGrid, subOrder, dimPoints(dimIdx), dimIdx, ~mod(dimIdx,2));
                    grid = vertcat(grid, newGrid); %#ok<AGROW>
                    order = vertcat(order, newOrder); %#ok<AGROW>
                end
                gridIdx = {gridIdx, order};
            end
            if nargout >= 2
                % Always select an xy slice, but make sure we map to the
                % correct order
                slice = xySlice(gridSize, gridIdx);
            end
        end
        
        % Diamond creator
        function [nPoints, slice, grid, order] = createDiamond(gridIdx, varargin)
            gridSize = meshGridSize(varargin{:});
            nPoints = ceil(prod(gridSize) / 2);
            if nargout >= (2 + iscell(gridIdx))
                % Start with a full mesh grid
                grids = cell(1, nargin-1);
                [grids{:}] = ndgrid(varargin{:});
                grid = flattenMerge(grids);
                order = columnColon(1,numel(grid));
                % Round dimensions up to next odd number.
                dimsOdd = 2.*floor(gridSize ./ 2) + 1;
                gridSel = false(dimsOdd);
                % Flood every other element.
                gridSel(1:2:end) = true;
                % Truncate dimensions that are even back to size.
                dimTrunc = arrayfun(@(d)colon(1,d),gridSize,'UniformOutput',false);
                gridSel = gridSel(dimTrunc{:});
                % Index into mesh grid to keep only required points.
                grid = grid(gridSel(:),:);
                % Diamond grids are in mesh grid order but with missing
                % entries.
                order = order(gridSel(:),1);
                gridIdx = {gridIdx, order};
            end
            if nargout >= 2
                % Always select an xy slice, but make sure we map to the
                % correct order
                slice = xySlice(gridSize, gridIdx);
            end
        end

        % Spiral creator
        function [nPoints, slice, grid, order] = createSpiralGrid(gridIdx, x, y, varargin)
            gridSize = meshGridSize(x, y, varargin{:});
            nPoints = prod(gridSize);
            if nargin == 3
                % Just X-Y, so simplify creation to a simple spiral.
                if nargout >= 3
                    grids = cell(1, 2);
                    [grids{:}] = ndgrid(x, y, varargin{:});
                    grid = flattenMerge(grids);
                    order = spiral(numel(x), numel(y));
                    grid = grid(order,:);
                end
                if nargout >= 2
                    % Slice is basically everything as already 2D
                    slice = defaultSlice(2);
                end
            else
                if nargout >= (2 + iscell(gridIdx))
                    % Get the sub-grid for all remaining dims
                    [~, ~, subGrid, subOrder] = GridPath.createSpiralGrid(gridIdx, x, y, varargin{1:end-1});
                    % Loop through all indicies in current dimension
                    dimPoints = varargin{end};
                    grid = [];
                    order = [];
                    for dimIdx = 1:numel(dimPoints)
                        % For each, add a new copy of the sub-grid for this new
                        % point. We flip the sub-grid for even points.
                        [newGrid, newOrder] = addNewDim(subGrid, subOrder, dimPoints(dimIdx), dimIdx, ~mod(dimIdx,2));
                        grid = vertcat(grid, newGrid); %#ok<AGROW>
                        order = vertcat(order, newOrder); %#ok<AGROW>
                    end
                    gridIdx = {gridIdx, order};
                end
                if nargout >= 2
                    % Always select an xy slice, but make sure we map to the
                    % correct order
                    slice = xySlice(gridSize, gridIdx);
                end
            end
        end
        
        % Spiral slice creator
        function [nPoints, slice, grid, order] = createSpiralSlice(gridIdx, x, y, z)
            import GridPath.createSpiralGrid
            narginchk(4,4);
            volSize = meshGridSize(x, y, z);
            nPoints = volSize(1)*volSize(2) + ...
                      volSize(1)*volSize(3) + ...
                      volSize(2)*volSize(3) - ...
                      sum(volSize) + 1;
            if nargout >= 3
                % Form each sub-spiral
                [nXY, ~, gridXY, orderXY] = createSpiralGrid(gridIdx, x, y);
                [nXZ, ~, gridXZ, orderXZ] = createSpiralGrid(gridIdx, x, z);
                [nYZ, ~, gridYZ, orderYZ] = createSpiralGrid(gridIdx, y, z);
                % Remap ordering.
                % - For XY it's simply shift the slice to the middle slice
                [orderXY_x,orderXY_y] = ind2sub(volSize([1 2]),orderXY);
                orderXY_zIdx = ceil(volSize(3)/2);
                orderXY_z = repmat(orderXY_zIdx, size(orderXY));
                orderXY = sub2ind(volSize, orderXY_x, orderXY_y, orderXY_z);
                gridXY = [gridXY(:,1), gridXY(:,2), repmat(z(orderXY_zIdx),nXY,1)];
                % - For XZ we need to shift to the middle slice and remap
                %   to the correct plane
                [orderXZ_x,orderXZ_z] = ind2sub(volSize([1 3]),orderXZ);
                orderXZ_yIdx = ceil(volSize(2)/2);
                orderXZ_y = repmat(orderXZ_yIdx, size(orderXZ));
                orderXZ = sub2ind(volSize, orderXZ_x, orderXZ_y, orderXZ_z);
                gridXZ = [gridXZ(:,1), repmat(y(orderXZ_yIdx),nXZ,1), gridXZ(:,2)];
                % - For YZ we need to shift to the middle slice and remap
                %   to the correct plane
                [orderYZ_y,orderYZ_z] = ind2sub(volSize([2 3]),orderYZ);
                orderYZ_xIdx = ceil(volSize(1)/2);
                orderYZ_x = repmat(orderYZ_xIdx, size(orderYZ));
                orderYZ = sub2ind(volSize, orderYZ_x, orderYZ_y, orderYZ_z);
                gridYZ = [repmat(x(orderYZ_xIdx),nYZ,1), gridYZ(:,1), gridYZ(:,2)];
                % Save the full grid and order
                grid = vertcat(gridXY, gridXZ, gridYZ);
                order = vertcat(orderXY, orderXZ, orderYZ);
                % Remove overlapping points from the cross.
                [grid, keep] = unique(grid,'rows','stable');
                order = order(keep);
            end
            if nargout >= 2
                % Determine which sector slice to display
                if iscell(gridIdx)
                    gridIdx = gridIdx{1};
                end
                boundaries = cumsum([volSize(1)*volSize(2), volSize(1)*volSize(3)-volSize(1)]);
                if isempty(gridIdx) || (gridIdx <= boundaries(1))
                    slice = {':', ':', ceil(volSize(3)/2)};
                elseif gridIdx <= boundaries(2)
                    slice = {':', ceil(volSize(2)/2), ':'};
                else
                    slice = {ceil(volSize(1)/2), ':', ':'};
                end
            end
        end

    end

end

% Ensure argument is correct size
function arg = checkArgSize(arg, nDims, logName)
if isscalar(arg)
    % Replicate if scalar
    arg = repmat(arg,1,nDims);
elseif numel(arg) < nDims
    % Else ensure large enough
    error('GridPath:ArgumentSize','Argument ''%s'' must be scalar or at least one element per dimension (>= %d)', logName, nDims);
else
    arg = reshape(arg(1:nDims),1,[]);
end

end

% Merge all points arrays, flattening them each first
function grid = flattenMerge(points)
points = cellfun(@(p)reshape(p,[],1), points, 'UniformOutput', false);
grid = horzcat(points{:});
end

% Determine full mesh grid size
function dims = meshGridSize(varargin)
dims = cellfun("prodofsize", varargin);
end

% Default slice array
function slice = defaultSlice(nDims)
slice = repmat({1},1,nDims);
[slice{1:min(2,end)}] = deal(':');
end

% X-Y slice array
function slice = xySlice(dimSize, gridIdx)
if iscell(gridIdx)
    gridIdx = gridIdx{2}(gridIdx{1});
end
slice = cell(1,numel(dimSize));
% Determine index in full grid for the specified index.
[slice{:}] = ind2sub(dimSize, gridIdx);
[slice{1:min(2,end)}] = deal(':');
end

% Add a new dimension to points array, with all values set to the same
% value, optionally flipping the original array first.
function [grid, order] = addNewDim(points, order, dimPoint, dimIdx, doFlip)
if doFlip
    grid = flipud(points);
    order = numel(order) * dimIdx + 1 - order;
else
    grid = points;
    order = order + numel(order) * (dimIdx - 1);
end
grid(:,end+1) = dimPoint;
end

% Colon as column
function col = columnColon(startVal,endOrStep,endVal)
if nargin < 3
    col = startVal:endOrStep;
else
    col = startVal:endOrStep:endVal;
end
col = col(:);
end

% NxM centre-out spiral
function [order] = spiral(nx, ny)

% First work corner in
nPts = nx * ny;
order = zeros(nPts,1);
sx = 1;sy = 2;
ex = nx;ey = ny;
curX = 1;
curY = 1;
dir = 0;
for idx = 1:nPts
    order(idx) = sub2ind([nx ny], curX, curY);
    switch dir
        case 0 % R
            if curX == ex
                dir = 1;
                ex = ex - 1;
                curY = curY + 1;
            else
                curX = curX + 1;
            end
        case 1 % D
            if curY == ey
                dir = 2;
                ey = ey - 1;
                curX = curX - 1;
            else
                curY = curY + 1;
            end
        case 2 % L
            if curX == sx
                dir = 3;
                sx = sx + 1;
                curY = curY - 1;
            else
                curX = curX - 1;
            end
        case 3 % U
            if curY == sy
                dir = 0;
                sy = sy + 1;
                curX = curX + 1;
            else
                curY = curY - 1;
            end
    end
end

% Reverse the order of the spiral to go outwards from centre
order = flipud(order);
end
