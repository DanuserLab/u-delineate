function varargout = FilamentAnalysisPackageGUI(varargin)
% Launch the GUI for the Filament Analysis Package
%
% This function calls the generic packageGUI function, passes all its input
% arguments and returns all output arguments of packageGUI
%
%
% Liya Ding  06/2012
%
% Updated by Qiongjing (Jenny) Zou, May 2024
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

if nargin>0 && isa(varargin{1},'MovieList')
    varargout{1} = packageGUI('FilamentAnalysisPackage',...
        [varargin{1}.getMovies{:}],'ML',varargin{1},varargin{2:end});
else
    varargout{1} = packageGUI('FilamentAnalysisPackage',varargin{:});
end
end