function tofData = filterInvalidInputs(tofData, amplitudeData,  amplitudeThreshold, tofMinValue, tofMaxValue, laplacianWindowSize, dx, dy)

% current assumption is that the grid is regular, the below would need to
% be updated if this was not the case:

if ~(isempty(dx) || isempty(dy) || isempty(laplacianWindowSize))
    for elementIdx = 1:size(tofData,3)
        subFilteredToFData = tofData(:, :, elementIdx);
        % Remove nan:
        subFilteredToFData(isnan(subFilteredToFData)) = 0;
        lap = del2(tofData(:, :, elementIdx), dx, dy);
        lapMag= abs(lap);
        expectedLap = medfilt2(lapMag, [laplacianWindowSize, laplacianWindowSize]);
        lapDiff = abs(lapMag - expectedLap);
        mask = lapDiff > 2 * expectedLap;
        subFilteredToFData(mask) = nan;
        tofData(:, :, elementIdx) = subFilteredToFData;
    end
end

if ~isempty(amplitudeData) && ~isempty(amplitudeThreshold)
    tofData(amplitudeData < amplitudeThreshold) = NaN;
end

if ~isempty(tofMinValue)
    tofData(tofData < tofMinValue) = NaN;
end

if ~isempty(tofMaxValue)
    tofData(tofData > tofMaxValue) = NaN;
end

end
