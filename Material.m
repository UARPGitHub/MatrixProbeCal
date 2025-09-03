classdef Material < handle
    % Base class for material helpers in UBase.Constant
    
    properties (Hidden, SetAccess = protected)
        SpeedOfSoundCoeff = zeros(1,7);
        DensityCoeff = zeros(1,7);
        TemperatureRange = [0 100];
    end
    
    methods (Hidden, Access = protected)
        % Constructor
        function obj = Material(sosCoeff, rhoCoeff, tempRange)
            obj.SpeedOfSoundCoeff = sosCoeff;
            obj.DensityCoeff = rhoCoeff;
            obj.TemperatureRange = tempRange;
        end
    end
    
    methods
        % Calculate speed of sound for fluid at temperature
        % This is an approximation over the specified valid temperature range.
        function [speedOfSound, density] = speedOfSound(obj, temperature)
            narginchk(2,2)
            
            minTemp = obj.TemperatureRange(1);
            maxTemp = obj.TemperatureRange(2);
            if any(temperature < minTemp) || any(temperature > maxTemp)
                error('Material:InvalidParameter','Approximations are only valid for range %d to %d %cC',minTemp, maxTemp, char(176));
            end
            
            % Calculate the speed of sound:
            speedOfSound = poly5calc(obj.SpeedOfSoundCoeff, temperature);
            
            % Calculate density
            density = poly5calc(obj.DensityCoeff, temperature);

        end
        
        % Calculate intensity from pressure
        % At a given fluid temperature, calculate the intensity from a peak
        % pressure measurement.
        % This is an approximation for 0*C to 90*C.
        function [intensity, speedOfSound, density] = pressureToIntensity(obj, pressure, temperature)
            narginchk(3,3)
            
            % Calculate the speed of sound and density
            [speedOfSound, density] = obj.speedOfSound(temperature);
            
            % Convert pressure to intensity
            intensity = (pressure.^2) ./ (density .* speedOfSound);
        end
    end
end

function res = poly5calc(coeffs, in)
% Calculates a nested fifth order polynomial with optional x^-1 term.
% Coefficients is a 7 term array with (1) being deepest nested coefficient,
% through to (6) being the constant term. (7) if non-zero is the reciprocal
% term.
res = (((((coeffs(1).*in + coeffs(2)).*in + coeffs(3)).*in + coeffs(4)).*in + coeffs(5)).*in + coeffs(6)) ./ (1 + coeffs(7).*in);
end

