% [res, theta, nms, scaleMap] = multiscaleSteerableDetector(img, M, varargin) performs multi-scale steerable filtering
% 
% Multi-scale filtering is performed by computing the steerable filter response at individual scales,
% and selecting the per-pixel maximum response across scales.
% The resulting scale map is used to produce an orientation map.
%
% Inputs: 
%            img : input image
%              M : order of the filter, between 1 and 5
%                : Odd orders: edge detectors, M = 1 is equivalent to Canny's detector
%                : Even orders: ridge detectors
%                  Higher orders provide better orientation selectivity and are less sensitive to noise,
%                  at a small trade-off in computational cost.
% {'sigmaArray'} : array of scales (standard deviation of the Gaussian kernel)
%                  across which the filters are evaluated. Default: [1 2 4]
%    
%
% Outputs: 
%         res : response to the filter across 
%       theta : orientation map
%         nms : non-maximum-suppressed response
%    scaleMap : pixel map of scales for which the response is maximal
%
% Recommended usage: 
%                    3rd/5th order for edge detection (M = 3 or 5)  
%                    4th order for ridge detection (M = 4)
%
% See also steerableDetector.m
%
% Copyright (C) 2025, Danuser Lab - UTSouthwestern 
%
% This file is part of FilamentAnalysisPackage.
% 
% FilamentAnalysisPackage is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% FilamentAnalysisPackage is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with FilamentAnalysisPackage.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

% Francois Aguet, 02/01/2012

function [res, theta, nms, scaleMap] = multiscaleSteerableDetector(img, M, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('img');
ip.addRequired('M', @(x) ismember(x, 1:5));
ip.addOptional('sigmaArray', [1 2 4]);
ip.parse(img, M, varargin{:});

sigma = ip.Results.sigmaArray;
ns = numel(sigma);

% arrays to store results from individual scales
ires = cell(1,ns);
itheta = cell(1,ns);
inms = cell(1,ns);

% responses at individual scales
for si = 1:ns
    [ires{si}, itheta{si}, inms{si}] = steerableDetector(img, ip.Results.M, sigma(si));
    % scale normalization for even orders
    if ~mod(M,2)
        ires{si} = sigma(si) * ires{si};
    end
end

% determine  f(x,y) = argmax_s r_s(x,y) (init. at scale 1)
res = ires{1};
theta = itheta{1};
scaleMap = ones(size(img));
for si = 2:ns
    idx = ires{si}>res;
    res(idx) = ires{si}(idx);
    theta(idx) = itheta{si}(idx);
    scaleMap(idx) = si;
end
nms = nonMaximumSuppression(res, theta);

