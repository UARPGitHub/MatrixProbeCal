function [elementCoordinates] = calculateElementLocations(cal, cnc, xdr)

% Get sorted CNC position order
cncPositions = cnc.Positions;
cncPositions(cnc.PositionOrder,:) = cncPositions;
cncPositions = reshape(cncPositions, [cnc.PositionGridSize, 3]);

% Calculate normalised measured amplitude
normalisedAmp = permute(cal.MeasuredAmplitude,[2 1]);
normalisedAmp = reshape(normalisedAmp ./ mean(normalisedAmp,'all','omitnan'), [cnc.PositionGridSize xdr.NTotalElements]);

% Filter measured distances to remove outliers
measuredDistance = reshape(permute(cal.MeasuredDistance,[2 1]), [cnc.PositionGridSize xdr.NTotalElements]);
measuredDistance = filterInvalidInputs(measuredDistance, normalisedAmp, ...
    cal.OutlierAmplitudeLimit, ...
    cal.OutlierDistanceLimit(1:min(end,1)), ...
    cal.OutlierDistanceLimit(2:min(end,2)), ...
    3, cnc.PositionGridStep(1), cnc.PositionGridStep(2));

% Experiment assumes transducer is nominally 10mm from
% hydrophone.
nominalZ = 10e-3;
elemPos = zeros(3, xdr.NTotalElements);

% Find best fit parabolas for each element
fprintf('Calculating raw element coordinates using LSQNONLIN. This may take a while...\n');
try
    lsqiter = 1e4;
    lsqtol = 1e-12;
    for elementIdx = 1:xdr.NTotalElements
        waitbar.next(1);
        % Filter first to ensure the function is differentiable:
        [elemData, cncData] = removeNaN(measuredDistance(:,:,elementIdx), cncPositions);
        if numel(elemData) < xdr.OutlierMinObservations
            % Assume dead element
            elemPos(:, elementIdx) = nan;
        else
            % Perform parabola best fit
            fun = @(xyz)euclidDistance(xyz, cncData, elemData);
            elemPos(:, elementIdx) = lsqnonlin(fun, ...
                        [median(cnc.PositionX), median(cnc.PositionY), nominalZ       ], ...
                        [   min(cnc.PositionX),    min(cnc.PositionY), nominalZ - 5e-3], ...
                        [   max(cnc.PositionX),    max(cnc.PositionY), nominalZ + 5e-3], ...
                        optimoptions('LSQNONLIN', ...
                                     'Display','off', ...
                                     'FunctionTolerance',lsqtol, ...
                                     'OptimalityTolerance',lsqtol, ...
                                     'StepTolerance',lsqtol, ...
                                     'MaxIterations', lsqiter ...
                                    ) ...
            );
        end
        if ~mod(elementIdx,10)
            fprintf('.');
        end
    end
    elemPos = reshape(elemPos,[3, xdr.NElements]);
catch
    warning('Failed to calculate element coordinates.');
end
fprintf('\n')
delete(waitbarCleanup);

% Plot the raw calculated locations to aid debugging
fig = figure();
clf(fig);
ax = axes(fig);
scatter3(ax, 1e3.*reshape(elemPos(1,:,:),[],1), ...
             1e3.*reshape(elemPos(2,:,:),[],1), ...
             1e3.*reshape(elemPos(3,:,:),[],1), ...
             'DisplayName','Raw Coordinates');
title(ax, 'Raw Element Coordinates');
xlabel(ax, 'Lateral (mm)');
ylabel(ax, 'Elevation (mm)');
zlabel(ax, 'Depth (mm)');

% Determine global panel rotation using rigid SVD transform
xdrIdealLateral = xdr.IdealElementLateralCoordinates;
xdrIdealElevation = xdr.IdealElementElevationCoordinates;
xdrIdealDepth = xdr.IdealElementDepthCoordinates;
xdrIdealLocations = permute(cat(3, xdrIdealLateral, xdrIdealElevation, xdrIdealDepth),[3 1 2]);

[globalRotationMatrix, globalTranslationMatrix] = rigidTransformSVD(elemPos, xdrIdealLocations);
rawElementPosn = reshape(elemPos, 3, []);
globalElementPosn = globalRotationMatrix.' * (rawElementPosn - globalTranslationMatrix);

fig = figure();
clf(fig);
ax = axes(fig);
scatter3(ax, 1e3.*(rawElementPosn(1, :) - globalTranslationMatrix(1)), ...
             1e3.*(rawElementPosn(2, :) - globalTranslationMatrix(2)), ...
             1e3.*(rawElementPosn(3, :) - globalTranslationMatrix(3)), ...
             'DisplayName','After Global Translation','MarkerEdgeColor','#77AC30');
title(ax, 'Element Coordinates');
legend(ax, 'Location', 'north');
xlabel(ax, 'Lateral (mm)');
ylabel(ax, 'Elevation (mm)');
zlabel(ax, 'Depth (mm)');
hold(ax,'on');
scatter3(ax, 1e3.*globalElementPosn(1, :), ...
             1e3.*globalElementPosn(2, :), ...
             1e3.*globalElementPosn(3, :), ...
             'DisplayName','After Global Rotation','MarkerEdgeColor','#D95319');

globalElementPosn = reshape(globalElementPosn, [3, xdr.NElements]);

% Determine sub-panel rotation
latSubPanelSize =   floor(xdr.NElements(1) ./ xdr.NSubPanels(1));
elevSubPanelSize = floor(xdr.NElements(2) ./ xdr.NSubPanels(2));
elementCoordinates = zeros(size(elemPos));
for elevSubIdx = 1:xdr.NSubPanels(2)
    % For each elevation sub-panel, process each lateral
    % sub-panel.
    for latSubIdx = 1:xdr.NSubPanels(1)
        elevRange = (1:elevSubPanelSize) + (elevSubIdx - 1) * elevSubPanelSize;
        latRange =  (1:latSubPanelSize ) + (latSubIdx  - 1) * latSubPanelSize ;
        xdrSubpanelIdealLocations = xdrIdealLocations(:, latRange, elevRange);
        subpanelElementPosn = globalElementPosn(:, latRange, elevRange);
        rawSubpanelElementPosn = reshape(xdrSubpanelIdealLocations, 3, []);
        
        [subpanelRotationMatrix, subpanelTranslationMatrix] = rigidTransformSVD(subpanelElementPosn, xdrSubpanelIdealLocations);
        subpanelElementPosn = (subpanelRotationMatrix * rawSubpanelElementPosn) + subpanelTranslationMatrix;
        
        elementCoordinates(:, latRange, elevRange) = reshape(subpanelElementPosn,[3 latSubPanelSize, elevSubPanelSize]);
    end
end
scatter3(ax, 1e3.*reshape(elementCoordinates(1, :),[],1), ...
             1e3.*reshape(elementCoordinates(2, :),[],1), ...
             1e3.*reshape(elementCoordinates(3, :),[],1), ...
             'DisplayName','After Panel Optimisation','MarkerEdgeColor','#A2142F');

elementCoordinates = reshape(elementCoordinates,[3, xdr.NElements]);

end
