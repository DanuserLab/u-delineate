function out_img = keep_largest_area(in_img)

labelMask = bwlabeln(in_img);

%Get their area
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
obAreas = regionprops(labelMask,'Area');
out_img=in_img;

%First, check that there are objects to remove
if length(obAreas) > 1
    obAreas = [obAreas.Area];
    %Sort by area
    [dummy,iSort] = sort(obAreas,'descend');
    %Keep only the largest requested number
    out_img = labelMask == iSort(1);
else
    out_img = in_img;
end
        