function combineChannelCellMaskCell = combineChannelMarkedCellAreaMask(MD)
% function to get cell mask from marked single cells, combining all channels
% Liya Ding, Jan, 2015
%
% Input:
%   MD:     The movieList object loaded before running this function
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

% initialize output, as a cell of all empty []
combineChannelCellMaskCell  = cell(1, MD.nFrames_);
for iFrame = 1 : MD.nFrames_
    combineChannelCellMaskCell{iFrame}=[];
end

% get the per channel per frame result
cellMaskCell = markedCellAreaMask(MD);


for iChannel = 1 : numel(MD.channels_)
    for iFrame = 1 : MD.nFrames_
        if(~isempty(cellMaskCell{iChannel, iFrame}))
            
            %if never anything, put zero image
            if(isempty(combineChannelCellMaskCell{iFrame}))
                combineChannelCellMaskCell{iFrame} = zeros(MD.imSize_);
            end
            
            % combine
            combineChannelCellMaskCell{iFrame} = (combineChannelCellMaskCell{iFrame} + cellMaskCell{iChannel, iFrame}) >0;
            
        end
    end
end
