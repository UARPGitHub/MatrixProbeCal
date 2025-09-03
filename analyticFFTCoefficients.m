function coeffs = analyticFFTCoefficients(coeffs, dim)
% Convert FFT coefficients to extract analytic signal

narginchk(1,2);

% Nothing to do if coefficients are empty.
if isempty(coeffs)
    return
end

% Determine size of FFT
if nargin < 2
    [~,dim] = max(size(coeffs));
end
n = size(coeffs, dim);
scaleSize = ones(1,ndims(coeffs));
scaleSize(dim) = n;

% Calculate scaling factors to extract analytic signal
scaling = zeros(scaleSize,'like',coeffs);
scaling([1 floor(n/2)+1]) = 1;
scaling(2:ceil(n/2)) = 2;

% Apply scaling
coeffs = scaling .* coeffs;

end
