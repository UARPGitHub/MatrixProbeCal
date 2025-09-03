function [rotationMat, translationMat] = rigidTransformSVD(rawPositions, idealPositions)

% Ensure flat [3 x Elem] arrays if multi-dimensional
rawPositions = reshape(rawPositions, 3, []);
idealPositions = reshape(idealPositions, 3, []);

% Remove NaNs:
deadElems = any(isnan(rawPositions),1);
idealPositions(:,deadElems) = [];
rawPositions(:,deadElems) = [];

% This function assumes the input arrays are (3 x N) where N is the
% number of elements:
centroidDetectedPosition = mean(rawPositions, 2);
centroidAssumedPosition = mean(idealPositions, 2);

% We need to translate both position matrices to centre around (0, 0, 0):
centeredDetectedPosition = rawPositions - centroidDetectedPosition;
centeredAssumedPosition = idealPositions - centroidAssumedPosition;

% Find the covariance matrix for the assumed and detected positions:
covMatrix = centeredAssumedPosition * (centeredDetectedPosition');

% Perform SVD to  compute the orthonormal matrices U and V:
[U, ~, V] = svd(covMatrix);

% Calculate the rotation matrix:
rotationMat = V * U';

if det(rotationMat) < 0
    % Recalulate the rotation matrix to avoid  reflection of the elements:
    V(:, 3) = -V(:, 3);
    rotationMat = V * U';
end

translationMat = centroidDetectedPosition - rotationMat * centroidAssumedPosition;

end
