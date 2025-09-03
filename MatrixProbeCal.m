%% Configure Parameters
% All units are SI, except temperatures which are in degrees Celcius.

% Grid Parameters
gridStep = 0.6e-3;
xSteps = 38;
ySteps = 41;

% CNC Grid
cnc = struct();
cnc.PositionX = GridPath.createUniformPointsCentred(0, gridStep, xSteps);
cnc.PositionY = GridPath.createUniformPointsCentred(0, gridStep, ySteps);
cnc.PositionGridSize = [xSteps, ySteps];
cnc.PositionGridStep = [gridStep gridStep];
cnc.PositionMesh = GridPath.SpiralGrid;
[cnc.Positions, ~, cnc.PositionOrder] = GridPath.createGrid(cnc.PositionMesh, cnc.PositionX, cnc.PositionY);

% XDR Parameters
xdr = struct();
xdr.NElements = [32 32];
xdr.NTotalElements = prod(xdr.NElements);
xdr.ElementPitch = [0.3e-3, 0.3e-3];
xdr.NSubPanels = [1, 4];
%  - Calculate element coordinates
segmentOffsets = symmetricArray(1,4,false) .* xdr.ElementPitch(2);
segmentOffsets = repmat(segmentOffsets,prod(xdr.NElements)/4,1);
[ elementIdxs ] = symmetricArray( 1, xdr.NElements, true );
[ elemLatGrid, elemElevGrid ] = ndgrid(elementIdxs{1} .* elementPitch(1), elementIdxs{2} .* elementPitch(2));
xdr.IdealElementDepthCoordinates = zeros(xdr.NElements);
xdr.IdealElementLateralCoordinates = elemLatGrid;
xdr.IdealElementElevationCoordinates = elemElevGrid + segmentOffsets;

% Calibration parameters
cal = struct();
%  - Test Medium
cal.Material = Fluid.Water;
%  - Minimum amplitude outlier filter
%    Will ignore any points with a recorded amplitide below this fraction peak amplitude
cal.OutlierAmplitudeLimit = 0.5;
%  - Bounds outlier filter
%    Will ignore any points with a distance outside this range. First
%    element is the lower bound, second is the upper bound. If set to
%    scalar value, will use lower bound only. If empty there will be no
%    limit set.
cal.OutlierDistanceLimit = 8e-3;
%  - Minimum number of data points to calculate an element position
%    If there are fewer points than this available after outliers are
%    removed, then the element will be ignored - this is useful for
%    avoiding dead elements affecting the panel best fit.
cal.OutlierMinObservations = 100;
%  - Some peak detection parameters to refine the peak detection:
%    - First is the level above the noise floor we want to filter the
%      peaks on:
cal.NoiseThresholdParam = 2;
%    - Second is the minimum prominence of the peak (prominence is based
%      on the height above a reference level):
cal.ProminenceThresholdParam = 0.05;
%    - Minimum peak width to be considered a peak as a fraction of an
%      ideal filtered and enveloped signal:
cal.PeakWidthThresholdParam = 0.3;
%  - Frequency resolution (the minimum delta for use in frequency domain)
cal.MinFrequencyResolution = 10e3;

% Oscilloscope Parameters
%  - Details of the oscilloscope configuration
scope = struct();
scope.NSamples = 32768;
scope.SampleRate = 800e6;

% Transmit Waveform
%  - Provide the transmit waveform that your system will excite each
%    element with to allow calculation of the matched filter.
%  - Ideally this should be a linear-frequency-modulated chirp signal to
%    allow for an accurate match.
transmit = struct();
transmit.SampleRate = 240e6;
transmit.Duration = 5e-6;
transmit.TimePoints = 0:1/transmit.SampleRate:transmit.Duration;
transmit.StartFrequency = 2.7e6;
transmit.StopFrequency = 4.5e6;
transmit.Waveform = chirp(transmit.TimePoints,transmit.StartFrequency, transmit.Duration, transmit.StopFrequency);

% Calculate the matched filter for the transmit waveform
%  - Determine FFT parameters for matched filtering
cal.FFTLength = 2^nextpow2(scope.SampleRate ./ cal.MinFrequencyResolution);
cal.FFTDelta = scope.SampleRate ./ cal.FFTLength;
% Resample the transmit waveform to match the receive sampling frequency:
filterWaveform = resample(transmitWaveform, scope.SampleRate, transmit.SampleRate);
% Perform the fft of transmit waveform:
fftFilterWaveform = squeeze(fft(filterWaveform, cal.FFTLength));
% Obtain the conjugate:
fftConjugateWaveform = conj(fftFilterWaveform);
% Calculate scaling factors to extract analytic signal
fftConjugateWaveform = analyticFFTCoefficients(fftConjugateWaveform);
cal.FFTFilterCoeff = reshape(fftConjugateWaveform, 1, 1, []);
% Find the 'ideal' width of the fft filtered pulse:
idealFilteredWaveform = abs(ifft(fftFilterWaveform .* fftConjugateWaveform, cal.FFTLength));
% For matched filter we need to shift the first element to the
% center of the array:
idealFilteredWaveform = fftshift(idealFilteredWaveform);
% Find the maximum peak width:
[~, location, width] = findpeaks(idealFilteredWaveform, ...
    'SortStr', "descend", ...
    'NPeaks', 1);
% To obtain the offset
cal.FilterOffset = location - floor(numel(idealFilteredWaveform)./2);
cal.PeakWidthThreshold = cal.PeakWidthThresholdParam .* width;


%% Collect Data
nPositions = numel(cnc.Positions);
cal.MeasuredAmplitude = nan(xdr.NTotalElements, nPositions);
cal.MeasuredDistance = nan(xdr.NTotalElements, nPositions);

% Loop through each position
for positionIdx = 1:nPositions
    
    % -----------------------------------------------------
    % Position CNC
    % Include your own code to move the CNC stage to the following position
    position = cnc.Positions(positionIdx);
    positionX = position(1);
    positionY = position(2);
    positionZ = position(3);
    
    % -----------------------------------------------------
    % Measure temperature of water
    % Insert own code for reading temperature from e.g. a thermocouple, or
    % use fixed value.
    temperature = 22; 
    
    % -----------------------------------------------------
    % Collect hydrophone data
    % Insert own code to fire each element of the transducer one by one.
    % For each element we capture hydrophone data using a scope. You can
    % use averaging (i.e. fire each element multiple times) to improve
    % hydrophone SNR.
    % This is all integrated into the UARPIII API, but any ultrasound setup
    % can be used.
    % For example:
    data = zeros(obj.NTotalElements, scope.NSamples);
    for elevIdx = 1:xdr.NElements(2)
        for latIdx = 1:xdr.NElements(1)
            % 1. Start scope capturing data
            % 2. Fire this element enough times for averages required
            % 3. Download averaged data from scope.
            data((elevIdx - 1) * xdr.NElements(1) + latIdx, :) = 0;
        end
    end
    
    % -----------------------------------------------------
    % Calculate the time of flight
    [distance, peak] = calculateTimeOfFlight(cal, xdr, scope, data);
    
    % -----------------------------------------------------
    % Save this data set to the correct place in point order
    pointIdx = cnc.PositionOrder(positionIdx);
    cal.MeasuredDistance(:,pointIdx) = distance;
    cal.MeasuredAmplitude(:,pointIdx) = peak;
    
end


%% Process the ToF Data

% Calculate the element locations
[elementCoordinates] = calculateElementLocations(cal, cnc, xdr);
xdr.ElementDepthCoordinates = elementCoordinates(1,:,:);
xdr.ElementLateralCoordinates = elementCoordinates(2,:,:);
xdr.ElementElevationCoordinates = elementCoordinates(3,:,:);



