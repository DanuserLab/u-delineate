function steerable_filter_forprocess_new(movieDataOrProcess, varargin)
% SteerableFilteringProcess wrapper function for SteerableFilteringProcess.
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
%
% This wrapper fcn was written based on steerable_filter_forprocess.m written by Liya Ding from 2015 or before.
% The way steerable_filter_forprocess was written is outdated, and not compatible with features on packageGUI, etc.
% Hillary Wong & Qiongjing (Jenny) Zou, Nov 2024
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
[movieData, thisProc] = getOwnerAndProcess(movieDataOrProcess, 'SteerableFilteringProcess', true);
p = parseProcessParams(thisProc, paramsIn); % If parameters are explicitly given, they should be used
% rather than the one stored in SteerableFilteringProcess

% Parameters: 

selected_channels = p.ChannelIndex;
BaseSteerableFilterSigma = p.BaseSteerableFilterSigma;
Levelsofsteerablefilters = p.Levelsofsteerablefilters;
ImageFlattenFlag = p.ImageFlattenFlag;
Sub_Sample_Num = p.Sub_Sample_Num;

nFrame = movieData.nFrames_;

SteerableFilteringProcessOutputDir = p.OutputDirectory;
currOutputDirectory = SteerableFilteringProcessOutputDir;

% Sanity Checks
nChan = numel(movieData.channels_);
if max(p.ChannelIndex) > nChan || min(p.ChannelIndex)<1 || ~isequal(round(p.ChannelIndex), p.ChannelIndex)
    error('Invalid channel numbers specified! Check ChannelIndex input!!')
end

% precondition / error checking
% Find the package of Filament Analysis
indexFilamentPackage = movieData.getPackageIndex('FilamentAnalysisPackage',1,true); % nDesired = 1 ; askUser = true
if isempty(indexFilamentPackage)
    error('Need to be in Filament Package for now.')
end

% Find ImageFlattenProcess
indexFlattenProcess = movieData.getProcessIndex('ImageFlattenProcess');
if isempty(indexFlattenProcess) && ImageFlattenFlag == 2
    error("The setting shows you want to use flattened image for steerable filtering. Please set parameters for Image Flatten and run.")
end

% logging input paths (bookkeeping)
% SteerableFilteringProcess uses the output of ImageFlattenProcess, if ImageFlattenFlag is 2, otherwise uses the raw images as input.
inFilePaths = cell(1, numel(movieData.channels_));
for i = p.ChannelIndex
    if ImageFlattenFlag ~= 2
        inFilePaths{1,i} = movieData.getChannelPaths{i};
    else
        inFilePaths{1,i} = movieData.processes_{indexFlattenProcess}.outFilePaths_{1,i};
    end
end
thisProc.setInFilePaths(inFilePaths);

% logging output paths.
dName = 'Channel';%String for naming the output directories for each channel
mkClrDir(currOutputDirectory, false);
outFilePaths = cell(1, numel(movieData.channels_));
for iChan = p.ChannelIndex
    % Create string for current directory
    currDir = [p.OutputDirectory filesep dName num2str(iChan)];
    outFilePaths{1,iChan} = currDir;
    mkClrDir(outFilePaths{1,iChan});
end
thisProc.setOutFilePaths(outFilePaths);

%% Algorithm
% see steerable_filter_forprocess.m line 117 and after

% RES_cell = cell(1,nFrame);
% Image_cell = cell(1,nFrame);

Frames_to_Seg = 1:Sub_Sample_Num:nFrame;
Frames_results_correspondence = im2col(repmat(Frames_to_Seg, [Sub_Sample_Num,1]),[1 1]);
Frames_results_correspondence = Frames_results_correspondence(1:nFrame);

for iChannel = selected_channels
    % Get frame number from the title of the image, this not neccesarily
    % the same as iFrame due to some shorting problem of the channel
    % Channel_FilesNames = movieData.channels_(iChannel).getImageFileNames(1:movieData.nFrames_);
    % QZ -- start
    fileNamesF = movieData.getImageFileNames(iChannel);
    Channel_FilesNames = fileNamesF{1}; % this is the way to get fileNames which works for diff kinds images, including tiffs, bioFormats images
    % QZ -- end

    filename_short_strs = uncommon_str_takeout(Channel_FilesNames);
    
    % Make output directory for the steerable filtered images
    % ImageSteerableFilterChannelOutputDir = movieData.processes_{indexSteerableProcess}.outFilePaths_{iChannel};
    ImageSteerableFilterChannelOutputDir = outFilePaths{1,iChannel}; % QZ

    % this part was done in logging output paths. - QZ
    % if (~exist(ImageSteerableFilterChannelOutputDir,'dir'))
    %     mkdir(ImageSteerableFilterChannelOutputDir);
    % end
    
    
    display('======================================');
    
    display(['Current movie: as in ',movieData.outputDirectory_]);
    
    display(['Start steerable filtering in Channel ',num2str(iChannel)]);
    
    if (~exist([ImageSteerableFilterChannelOutputDir,filesep,'NMS'],'dir'))
        try
            mkdir(ImageSteerableFilterChannelOutputDir,'NMS');
        catch
            system(['mkdir -p ' ImageSteerableFilterChannelOutputDir filesep 'NMS']);
        end
    end
    
    for iFrame_subsample = 1 : length(Frames_to_Seg)
        iFrame = Frames_to_Seg(iFrame_subsample);
        disp(['Frame: ',num2str(iFrame)]);

        TIC_IC_IF = tic;
        % Read in the intensity image.
        if indexFlattenProcess > 0
            currentImg = imread([movieData.processes_{indexFlattenProcess}.outFilePaths_{iChannel}, filesep, 'flatten_',filename_short_strs{iFrame},'.tif']);
        else
            currentImg = movieData.channels_(iChannel).loadImage(iFrame);
        end
        currentImg = double(currentImg);
        
        levels_sizes = 2.^((1:Levelsofsteerablefilters)-1);
        
        % Steerable filtering using four scales one doubling the previous one.
        % function multiscaleSteerableDetector will automatically merge the results
%         [MAX_st_res, orienation_map, nms, scaleMap] = multiscaleSteerableDetector(currentImg, 4, BaseSteerableFilterSigma.*levels_sizes);
        
        [MAX_st_res, orienation_map, nms, scaleMap] = multiscaleSteerableDetector(currentImg, 4, BaseSteerableFilterSigma.*levels_sizes);
        
        for sub_i = 1 : Sub_Sample_Num
            if iFrame + sub_i-1 <= nFrame
                imwrite((MAX_st_res)/(max(max(MAX_st_res))), ...
                    [ImageSteerableFilterChannelOutputDir,filesep,'MAX_st_res_', ...
                    filename_short_strs{iFrame + sub_i-1},'.tif']);
               imwrite((nms)/(max(max(nms))), ...
                    [ImageSteerableFilterChannelOutputDir,filesep,'NMS',filesep,'NMS_', ...
                    filename_short_strs{iFrame + sub_i-1},'.tif']);
            end
        end
        
        % change these into single to save space in drug screen
        % for most experiments, single precision is enough.
        orienation_map = single(orienation_map);
        MAX_st_res = single(MAX_st_res);
        nms = single(nms);
        scaleMap = single(scaleMap);
        
        save([ImageSteerableFilterChannelOutputDir,filesep,'steerable_',filename_short_strs{iFrame},'.mat'],...
            'orienation_map', 'MAX_st_res','nms','scaleMap');
        
        Time_cost = toc(TIC_IC_IF);
        disp(['Frame ', num2str(iFrame), ' ST filter costed ',num2str(Time_cost,'%.2f'),'s.']);                        

    end
end
%%%% end of algorithm

fprintf('Finished Steerable Filtering Process! \n')
end

