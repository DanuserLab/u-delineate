function output_strs = all_uncommon_str_takeout(input_strs)

% old fashion, when only one file, only the last character is kept.
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
if iscell(input_strs)==0
    output_strs = input_strs(end);
    return;
end


stop_num = common_str_stop(input_strs);
cutback_num = common_str_cutback(input_strs);



for i_str = 1 : length(input_strs)
    output_strs{i_str} = input_strs{i_str}(stop_num:end-cutback_num);
end
