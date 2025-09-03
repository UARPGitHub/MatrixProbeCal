classdef Fluid < Material
    % Fluid Properties
    
    enumeration
        % Different types of fluid along with coefficients for calculating
        % various parameters at different temperatures.
        % First argument is speed of sound coefficients for nested form 5th
        % order polynomial with x^-1 term.
        % Second argument is density coefficients in same polynomial format.
        %       [X5          ,X4          ,X3          ,X2          ,X1          ,X0          ,X-1         ]
        
        % Air at standard atmospheric pressure
        Air   ( [ 0          , 0          , 9.150466e-7,-5.972844e-4, 0.605386295, 331.423825 , 0          ], ...  % Speed of Sound (Engineering Toolbox)
                [ 0          , 0          ,-4.539475e-8, 1.760323e-5,-4.763127e-3, 1.29223502 , 0          ], ...  % Density (Engineering Toolbox)
                [-40 100]); % Valid temperature range
        
        % Water at standard atmospheric pressure
        Water ( [ 3.16585e-9 ,-1.4826e-6  , 3.34638e-4 ,-5.81173e-2 , 5.03836    , 1402.39    , 0          ], ...  % Speed of Sound (Bilaniuk and Wong, 148pt Equation),
                [-2.80543e-10, 1.05563e-7 ,-4.617046e-5,-7.98704e-3 , 16.945176  , 999.83952  , 1.687985e-2], ...  % Density (Kell et al.)
                [0 90]); % Valid temperature range
    end
    
    methods (Hidden, Access = protected)
        % Constructor
        function obj = Fluid(sosCoeff, rhoCoeff, tempRange)
            obj = obj@Material(sosCoeff, rhoCoeff, tempRange);
        end
    end
end
