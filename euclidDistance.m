function euclidXYZ = euclidDistance(xyz0, cncPos, elemDist)
% Calculate euclidean distance for lsqnonlin

euclidXYZ = sqrt((cncPos(:, 1)-xyz0(1)).^2 + (cncPos(:, 2)-xyz0(2)).^2 + (cncPos(:, 3)-xyz0(3)).^2) - elemDist;

end
