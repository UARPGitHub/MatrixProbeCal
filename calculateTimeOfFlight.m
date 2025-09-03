function [distance, peak] = calculateTimeOfFlight(cal, xdr, scope, data)

% Apply matched filtering and calculate time of flight
data = fft(data, cal.FFTLength, 2);
data = data .* cal.FFTFilterCoeff;
data = abs(ifft(data, cal.FFTLength, 2));
data = data(:, 1:scope.NSamples);

% Calculate the speed of sound for the test medium;:
speedOfSound = cal.Material.speedOfSound(temperature);

% Calculate the median of the input data as this can give a
% noise floor threshold.
noiseThreshold = cal.NoiseThresholdParam .* median(data, 2);

% Naive peak detection used to get the amplitude for prominence
% thresholding:
[naivePeak, ~] = max(data, [], 2);
prominenceThreshold = cal.ProminenceThresholdParam .* naivePeak;

% Calculate the peak, unfortunately, the matlab peak detection
% only works on a vector so we need a loop here:
peak = nan(xdr.NTotalElements,1);
peakIndex = ones(xdr.NTotalElements,1);
warnStruct = warning('off','signal:findpeaks:largeMinPeakHeight');
for elementItr = 1:xdr.NTotalElements

    % Only ever return one peak
    [peakValue, peakIndexValue] = findpeaks(data(elementItr, :), ...
        'MinPeakHeight', noiseThreshold(elementItr), ...
        'MinPeakProminence', prominenceThreshold(elementItr), ...
        'MinPeakWidth', cal.PeakWidthThreshold, ...
        'SortStr', "descend", ...
        'NPeaks', 1);

    % If we managed to find a peak based on our thresholding
    % attempt then store the peak for later distance
    % calculations:
    if ~isempty(peakIndexValue)
        peak(elementItr) = peakValue;
        peakIndex(elementItr) = peakIndexValue;
    end
end
warning(warnStruct);

% And convert peak index to position:
distance = speedOfSound * ((peakIndex-cal.FilterOffset)/scope.SampleRate);

end
