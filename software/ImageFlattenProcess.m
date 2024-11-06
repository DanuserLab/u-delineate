classdef ImageFlattenProcess < ImageProcessingProcess
    % A concrete class for flatten images
    %
    % Liya Ding, 06. 2012
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
    
    methods (Access = public)
        
        function obj = ImageFlattenProcess(owner,varargin)            
            
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = ImageFlattenProcess.getName;
                super_args{3} = @image_flatten_new; % edited 2024-10. The wrapper fcn was image_flatten.
                if isempty(funParams)
                    funParams = ImageFlattenProcess.getDefaultParams(owner,outputDir); % Added outputDir here so the result folder is saved in the package's folder. - QZ Nov 2024
                end
                super_args{4} = funParams;
                
                if nargin > 4
                    super_args{5} = inImagePaths;
                end
                if nargin > 5
                    super_args{6} = outImagePaths;
                end
            end
            
            obj = obj@ImageProcessingProcess(super_args{:});
            
        end
                
        
         function setInImagePath(obj,chanNum,imagePath)
            
            if ~obj.checkChanNum(chanNum)
                error('lccb:set:fatal','Invalid image channel number for image path!\n\n');
            end
            
            if ~iscell(imagePath)
                imagePath = {imagePath};
            end
            nChan = length(chanNum);
            if nChan ~= length(imagePath)
                error('lccb:set:fatal','You must specify a path for every channel!')
            end
            
            for j = 1:nChan
                if ~exist(imagePath{j},'dir')
                    error('lccb:set:fatal',...
                        ['The directory specified for channel ' ...
                        num2str(chanNum(j)) ' is invalid!'])
                    
                else
                    if isempty(imDir(imagePath{j})) && ...
                            isempty(dir([imagePath{j} filesep '*.mat']))
                        error('lccb:set:fatal',...
                            ['The directory specified for channel ' ...
                            num2str(chanNum(j)) ' does not contain any images!!'])
                    else
                        obj.inFilePaths_{1,chanNum(j)} = imagePath{j};
                    end
                end
            end
         end
         
         function fileNames = getOutImageFileNames(obj,iChan)
            if obj.checkChannelOutput(iChan)
                fileNames = cellfun(@(x)(dir([x filesep '*.tif'])),obj.outFilePaths_(1,iChan),'UniformOutput',false);
                fileNames = cellfun(@(x)(arrayfun(@(x)(x.name),x,'UniformOutput',false)),fileNames,'UniformOutput',false);
                nChan = numel(iChan);
                for j = 1:nChan
                    %Sort the files by the trailing numbers
                    fNums = cellfun(@(x)(str2double(...
                        x(max(regexp(x(1:end-4),'\D'))+1:end-4))),fileNames{j});
                    [~,iX] = sort(fNums);
                    fileNames{j} = fileNames{j}(iX);
                end
                nIm = cellfun(@(x)(length(x)),fileNames);
                if ~all(nIm == obj.owner_.nFrames_)
                    error('Incorrect number of images found in one or more channels!')
                end
            else
                error('Invalid channel numbers! Must be positive integers less than the number of image channels!')
            end
            
            
         end
        
        
            function fileNames = getInImageFileNames(obj,iChan)
            if obj.checkChanNum(iChan)
                
                nChan = numel(iChan);
                fileNames = cell(1,nChan);
                for j = 1:nChan
                    %First check for regular image inputs
                    fileNames{j} = imDir(obj.inFilePaths_{1,iChan(j)});
                    if isempty(fileNames{j})
                        %If none found, check for .mat image inputs
                        fileNames{j} = dir([obj.inFilePaths_{1,inFilePaths_iChan(j)} filesep '*.tif']);
                    end
                    fileNames{j} = arrayfun(@(x)(x.name),fileNames{j},'UniformOutput',false);
                    nIm = length(fileNames{j});
                    if nIm ~= obj.owner_.nFrames_
                        error(['Incorrect number of images found in channel ' num2str(iChan(j)) ' !'])
                    end
                end
            else
                error('Invalid channel numbers! Must be positive integers less than the number of image channels!')
            end
            
            
            end
        
            
            
        function setOutImagePath(obj,chanNum,imagePath)
            
            if ~obj.checkChanNum(chanNum)
                error('lccb:set:fatal','Invalid image channel number for image path!\n\n');
            end
            
            if ~iscell(imagePath)
                imagePath = {imagePath};
            end
            nChan = length(chanNum);
            if nChan ~= length(imagePath)
                error('lccb:set:fatal','You must specify a path for every channel!')
            end
            
            for j = 1:nChan
                if ~exist(imagePath{j},'dir')
                    error('lccb:set:fatal',...
                        ['The directory specified for channel ' ...
                        num2str(chanNum(j)) ' is invalid!'])
                else
                    obj.outFilePaths_{1,chanNum(j)} = imagePath{j};
                end
            end
            
            
        end
        
        
        
         function out_data = loadChannelOutput(obj,iChan,iFrame,varargin)
            % Input check
            ip =inputParser;
            ip.addRequired('iChan',@obj.checkChanNum);
            ip.addRequired('iFrame',@obj.checkFrameNum);
            
            outputList = {'flattened_image',''};
            ip.addParamValue('output',{},@(x) all(ismember(x,outputList)));
            
            ip.parse(iChan,iFrame,varargin{:})
            
            ImageFlattenChannelOutputDir = obj.outFilePaths_{iChan};
    
            Channel_FilesNames = obj.getOutImageFileNames(iChan); % Use output images file names instead of input image file names, to fixed the bug if  BioFormats data as input. - Qiongjing (Jenny) Zou, Nov 2024
            filename_short_strs = uncommon_str_takeout(Channel_FilesNames{1});
            
            % this line in commandation for shortest version of filename
            filename_shortshort_strs = all_uncommon_str_takeout(Channel_FilesNames{1});
            
            currentImg=[];
            
            try
                currentImg = imread([ImageFlattenChannelOutputDir,filesep, ... 
                    filename_short_strs{iFrame},'.tif']);
            catch
                try
                currentImg = imread([ImageFlattenChannelOutputDir,filesep, ...
                    filename_shortshort_strs{iFrame},'.tif']);
                catch
                    if(iFrame==1)
                        try
                            currentImg = imread([ImageFlattenChannelOutputDir,filesep,'flatten_f.tif']);
                        catch
                            currentImg = imread([ImageFlattenChannelOutputDir,filesep,'flatten_F.tif']);
                        end
                    end
                end
            end
            
            out_data = currentImg;
            
         end         
        
        function h = draw(obj,iChan,varargin)
            
            outputList = obj.getDrawableOutput();
            drawImageFlattenImage = any(strcmpi('ImageFlatten',varargin));
            
            if drawImageFlattenImage
                % Input check
                ip =inputParser;
                ip.addRequired('iChan',@(x) ismember(x,1:numel(obj.owner_.channels_)));
                ip.addParamValue('output',[],@ischar);
                ip.KeepUnmatched = true;
                ip.parse(iChan,varargin{:})
                
                % Load average corrected image
                s = load(obj.outFilePaths_{2,iChan});
                tmpFields=fieldnames(s);
                data=s.(tmpFields{1});
                
                iOutput= find(cellfun(@(y) isequal(ip.Results.output,y),{outputList.var}));
                if ~isempty(outputList(iOutput).formatData),
                    data=outputList(iOutput).formatData(data);
                end
                
                try
                    assert(~isempty(obj.displayMethod_{iOutput,iChan}));
                catch ME
                    obj.displayMethod_{iOutput,iChan}=...
                        outputList(iOutput).defaultDisplayMethod(iChan);
                end
                
                % Delegate to the corresponding method
                tag = [obj.getName '_channel' num2str(iChan) '_output' num2str(iOutput)];
                drawArgs=reshape([fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]',...
                    2*numel(fieldnames(ip.Unmatched)),1);
                h=obj.displayMethod_{iOutput,iChan}.draw(data,tag,drawArgs{:});
            else
                h=draw@ImageProcessingProcess(obj,iChan,varargin{1},varargin{2:end});
            end
        end
    end

    
    methods (Static)
        function name =getName()
            name = 'Image Flatten';
        end
        function h = GUI()
            h= @imageFlattenProcessGUI;
        end
        
        function output = getDrawableOutput()
            output = ImageProcessingProcess.getDrawableOutput();
        end
        
        function methods = getMethods(varargin)
            flatteningMethods(1).name = 'Log';
            flatteningMethods(2).name = 'Sqrt';
            flatteningMethods(3).name = 'Power 2/3';
            
            ip=inputParser;
            ip.addOptional('index',1:length(flatteningMethods),@isvector);
            ip.parse(varargin{:});
            index = ip.Results.index;
            methods=flatteningMethods(index);
        end
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_, @ischar);
            ip.parse(owner, varargin{:})
            outputDir = ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1:numel(owner.channels_);
            funParams.OutputDirectory = [outputDir  filesep 'ImageFlatten']; % Added 2024-10 to fix folder icon issue on packageGUI
%             funParams.outputDir = outputDir;
            funParams.method_ind = 3;
            funParams.imageflattening_mode = 2;
            funParams.GaussFilterSigma = 0.2;                        
            funParams.TimeFilterSigma = 0;    
            
            funParams.stat.low_005_percentile = 0;
            funParams.stat.high_995_percentile = 2^16-1;
            funParams.stat.center_value_int = 300;

            % sub-sample number, since often VIF images are taken at a
            % lower sample rate than the other channel, so use this number
            % to save some time.
            funParams.Sub_Sample_Num = 1;

        end
    end
end