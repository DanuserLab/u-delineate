function stop_num = common_str_stop(input_strs)

if length(input_strs)<=1
    stop_num =1;
    return;
end

if iscell(input_strs)==0
    stop_num =1;
    return;
end

reached = 0;

for stop_num = 1 : length(input_strs{1})    
    for i_str = 2 : length(input_strs)
        try
            if(~strcmp(input_strs{1}(1:stop_num),input_strs{i_str}(1:stop_num)))
                reached=1;
                break;
            end
        catch
            reached=1;
            break;
        end
    end
    if(reached==1)
        break;
    end
end
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
