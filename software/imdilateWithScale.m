function bw_out = imdilateWithScale(bw_in,scaleMap,levels)
% dilate a bw image with the corresponding scale
%
% Copyright (C) 2024, Danuser Lab - UTSouthwestern 
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

bw_out = bw_in;

for i = 1 : length(levels)
    if levels(i)>0
        H = fspecial('disk', ceil(levels(i)))>0;
        this_scale_bw = and(scaleMap==i, bw_in>0);
        this_scale_bw = imdilate(this_scale_bw, H,'same');
        bw_out = or(bw_out,this_scale_bw);
    end
end
