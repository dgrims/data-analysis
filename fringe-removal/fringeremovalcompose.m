function [optrefimages] = fringeremovalcompose(refimages,optcoeff)
% fringeremovalcompose creates an optimal reference image for each absorption image
% in a set as a linear combination of reference images, using the supplied
% coefficients. These coefficents have to be calculated using fringeremoval.m
%
% Input:
%   refimages = Raw reference image data
%   optcoeff = Array of optimal coefficents for each image.
%
% Output:
%   optrefimages = Array of optimal reference images
%
%   Version 1.0.0   Vladislav Gavryusev
%   History:
%       v1.0.0 26/09/2016 First version.

% Set variables, and flatten reference images
[ydim,xdim,nimgsR] = size(refimages);
nimgs = size(optcoeff,2);

R = single(reshape(refimages,xdim*ydim,nimgsR));
optrefimages = zeros([ydim,xdim,nimgs]); % preallocate

for j=1:nimgs % Compute optimised reference image
    optrefimages(:,:,j) = reshape(R*optcoeff(:,j),[ydim xdim]);
end

end