function RGB_seg_orient_heat_map = heat_display_filament_from_model(currentImg, current_model)

[digital_model,orientation_model,VIF_XX, VIF_YY,VIF_OO]...
    = filament_model_to_digital_with_orientation(current_model);
currentImg(currentImg==0)=46;

currentImg = 255*(1-(currentImg-min(min(currentImg)))/(max(max(currentImg))-min(min(currentImg))));
currentImg = uint8(currentImg);
img_size = size(currentImg);

Hue = (VIF_OO(:)+pi/2)/(pi)-0.2;
Hue(find(Hue>=1)) = Hue(find(Hue>=1)) -1;
Hue(find(Hue<0)) = Hue(find(Hue<0)) +1;

Sat = Hue*0+1;
Value = Hue*0+1;

RGB_seg_orient_heat_array = hsv2rgb([Hue Sat Value]);

R_seg_orient_heat_map = currentImg;
G_seg_orient_heat_map = currentImg;
B_seg_orient_heat_map = currentImg;

R_seg_orient_heat_map(sub2ind(img_size,VIF_YY,VIF_XX))=255*RGB_seg_orient_heat_array(:,1);
G_seg_orient_heat_map(sub2ind(img_size,VIF_YY,VIF_XX))=255*RGB_seg_orient_heat_array(:,2);
B_seg_orient_heat_map(sub2ind(img_size,VIF_YY,VIF_XX))=255*RGB_seg_orient_heat_array(:,3);

RGB_seg_orient_heat_map = uint8(zeros(img_size(1),img_size(2),3));

RGB_seg_orient_heat_map(:,:,1) = R_seg_orient_heat_map;
RGB_seg_orient_heat_map(:,:,2) = G_seg_orient_heat_map;
RGB_seg_orient_heat_map(:,:,3) = B_seg_orient_heat_map;

% figure;imshow(RGB_seg_orient_heat_map);
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

