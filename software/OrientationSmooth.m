function O_filter = OrientationSmooth(input_Orientation,SteerabelRes_Segment)

DO = input_Orientation*2;
DA = cos(DO);
DB = -sin(DO);
DA_filter = DA;
DB_filter = DB;

radius = 3;
H = fspecial('gaussian',2*radius+1,3);

for i = radius+1 : size(DA_filter,1)-radius
    for j = radius+1 : size(DA_filter,2)-radius
        if(SteerabelRes_Segment(i,j)==1)
        local_seg = SteerabelRes_Segment(i-radius:i+radius, j-radius:j+radius);
        local_H = local_seg.*H;
        local_H = local_H/(sum(local_H(:)));
        
        local_A_patch = DA(i-radius:i+radius, j-radius:j+radius);
        local_B_patch = DB(i-radius:i+radius, j-radius:j+radius);
%         local_A_sm = imfilter(local_A_patch,local_H,'replicate','same');
%         local_B_sm = imfilter(local_B_patch,local_H,'replicate','same');
%         
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
        DA_filter(i,j) = sum(sum(local_A_patch.*local_H));
        DB_filter(i,j) = sum(sum(local_B_patch.*local_H));

        TT = [DA_filter(i,j) DB_filter(i,j)];
        if abs(norm(TT)-1)>0.01
            TT = TT/norm(TT);
        end
        DA_filter(i,j) = TT(1);
        DB_filter(i,j) = TT(2);        

        
        end
    end
end

% DA_filter = imfilter(DA,fspecial('gaussian',11,5),'replicate','same'); 
% DB_filter = imfilter(DB,fspecial('gaussian',11,5),'replicate','same');
% 
% for i = 1 : size(DA_filter,1)
%     for j = 1 : size(DA_filter,2)
%         TT = [DA_filter(i,j) DB_filter(i,j)];
%         if abs(norm(TT)-1)>0.01
%             TT = TT/norm(TT);
%         end
%         DA_filter(i,j) = TT(1);
%         DB_filter(i,j) = TT(2);        
%     end
% end

DO_angle = atan2(-DB_filter,DA_filter);
O_filter = DO_angle/2;