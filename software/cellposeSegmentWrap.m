function cellposeSegmentWrap(movieDataOrProcess, varargin)
% cellposeSegmentWrap wrapper function for CellposeSegmentationProcess.
%
% INPUT
% movieDataOrProcess - either a MovieData (legacy)
%                      or a Process (new as of July 2016)
%
% param - (optional) A struct describing the parameters, overrides the
%                    parameters stored in the process (as of Aug 2016)
%
% OUTPUT
% none (saved to p.OutputDirectory)
%
% Changes
% As of July 2016, the first argument could also be a Process. Use
% getOwnerAndProcess to simplify compatability.
%
% As of August 2016, the standard second argument should be the parameter
% structure
%
% Qiongjing (Jenny) Zou, June 2024
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

%% ------------------ Input ---------------- %%
ip = inputParser;
ip.addRequired('MD', @(x) isa(x,'MovieData') || isa(x,'Process') && isa(x.getOwner(),'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieDataOrProcess, varargin{:});
paramsIn = ip.Results.paramsIn;

%% Registration
% Get MovieData object and Process
[movieData, thisProc] = getOwnerAndProcess(movieDataOrProcess, 'CellposeSegmentationProcess', true);
p = parseProcessParams(thisProc, paramsIn); % If parameters are explicitly given, they should be used
% rather than the one stored in CellposeSegmentationProcess

% Parameters:
% p

% Sanity Checks
nChan = numel(movieData.channels_);
if max(p.ChannelIndex) > nChan || min(p.ChannelIndex)<1 || ~isequal(round(p.ChannelIndex), p.ChannelIndex)
    error('Invalid channel numbers specified! Check ChannelIndex input!!')
end

% precondition / error checking
% check if ImagePreprocessingProcess was run
if isempty(p.ProcessIndex)
    iImagePreprocessProc = movieData.getProcessIndex('ImagePreprocessingProcess',1,true); % nDesired = 1 ; askUser = true
    if isempty(iImagePreprocessProc)
        error('ImagePreprocessingProcess needs to be done before run this process.')
    end
elseif isa(movieData.processes_{p.ProcessIndex},'ImagePreprocessingProcess')
    iImagePreprocessProc = p.ProcessIndex;
else
    error('The process specified by ProcessIndex is not a valid ImagePreprocessingProcess! Check input!')
end

% checking if GPU available, if not, give warning. GPU is highly recommended for speed!
if gpuDeviceCount("available") < 1
    fprintf(2, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
    warning('No GPU found! GPU is highly recommended for Cellpose Segmentation Process for speed!')
    fprintf(2, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
end

% logging input paths (bookkeeping)
inFilePaths = cell(1, numel(movieData.channels_));
for i = p.ChannelIndex
   inFilePaths{1,i} = movieData.processes_{iImagePreprocessProc}.outFilePaths_{1,i}; 
end
thisProc.setInFilePaths(inFilePaths);

% logging output paths.
mkClrDir(p.OutputDirectory);
outFilePaths = cell(1, numel(movieData.channels_));
for i = p.ChannelIndex
    outFilePaths{1,i} = [p.OutputDirectory filesep 'ch' num2str(i)];
    mkClrDir(outFilePaths{1,i});
end
thisProc.setOutFilePaths(outFilePaths);


%% Algorithm
% see segment3D_cellpose_segmentMD.py
pyenv(ExecutionMode="OutOfProcess")

% When pass Matlab variables as inputs in pyrunfile, it automatically converts MATLAB data into types that best represent the data to the Python language.
% Such as logical in Matlab to bool in Python, double in Matlab to float in Python.
% see more https://www.mathworks.com/help/matlab/matlab_external/passing-data-to-python.html

for k = p.ChannelIndex
    inPath = inFilePaths{1,k};
    outPath = outFilePaths{1,k};

    pPassToPython = rmfield(p,{'ChannelIndex','OutputDirectory','ProcessIndex'}); % only pass the same parameters in cellpose_segment_params to Python script

    pPassToPython.saveplotsfolder = outPath;

    pyrunfile("segment3D_cellpose_segmentMD.py", input_path = inPath, output_path = outPath, ...
        params = pPassToPython) % pyrunfile Since R2021b
end

%%%% end of algorithm

fprintf('\n Finished Cellpose Segmentation Process! \n')

end