function image_flatten_new(movieDataOrProcess, varargin)
% ImageFlattenProcess wrapper function for ImageFlattenProcess.
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
% This wrapper fcn was written based on image_flatten.m written by Liya Ding from 2015 or before.
% The way image_flatten was written is outdated, and not compatible with features on packageGUI, etc.
% Qiongjing (Jenny) Zou, Oct 2024
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
[movieData, thisProc] = getOwnerAndProcess(movieDataOrProcess, 'ImageFlattenProcess', true);
p = parseProcessParams(thisProc, paramsIn); % If parameters are explicitly given, they should be used
% rather than the one stored in ImageFlattenProcess

% Parameters: 
selected_channels = p.ChannelIndex;
flatten_method_ind = p.method_ind;
Gaussian_sigma = p.GaussFilterSigma;

TimeFilterSigma = p.TimeFilterSigma;
Sub_Sample_Num  = p.Sub_Sample_Num;

imageflattening_mode = p.imageflattening_mode;

nFrame = movieData.nFrames_;

ImageFlattenProcessOutputDir = p.OutputDirectory;
currOutputDirectory = ImageFlattenProcessOutputDir;

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

% logging input paths (bookkeeping)
% ImageFlattenProcess uses the raw images as input
inFilePaths = cell(1, numel(movieData.channels_));
for i = p.ChannelIndex
    inFilePaths{1,i} = movieData.getChannelPaths{i};
end
thisProc.setInFilePaths(inFilePaths);

% logging output paths.
dName = 'Channel';%String for naming the output directories for each channel
mkClrDir(currOutputDirectory, false);
outFilePaths = cell(1, numel(movieData.channels_));
for iChan = p.ChannelIndex;
    % Create string for current directory
    currDir = [p.OutputDirectory filesep dName num2str(iChan)];
    outFilePaths{1,iChan} = currDir;
    mkClrDir(outFilePaths{1,iChan});
end
thisProc.setOutFilePaths(outFilePaths);


%% Algorithm
% see image_flatten.m line 104 and after

% for most cases, don't do background removal here.
background_removal_flag =0;

for iChannel = selected_channels
    display('======================================');
    display(['Current movie: as in ',movieData.outputDirectory_]);

    for_temporal_filtering = cell(1);
    
    % ImageFlattenProcessOutputDir = movieData.processes_{indexFlattenProcess}.outFilePaths_{iChannel};
    ImageFlattenChannelOutputDir = outFilePaths{1,iChannel}; % QZ

    % Get frame number from the title of the image, this not neccesarily
    % the same as iFrame due to some shorting problem of the channel
    % Channel_FilesNames = movieData.channels_(iChannel).getImageFileNames(1:movieData.nFrames_);
    % QZ -- start
    fileNamesF = movieData.getImageFileNames(iChannel); 
    Channel_FilesNames = fileNamesF{1}; % this is the way to get fileNames which works for diff kinds images, including tiffs, bioFormats images
    % QZ -- end
    filename_short_strs = uncommon_str_takeout(Channel_FilesNames);
    
    Frames_to_Seg = 1:Sub_Sample_Num:nFrame;
    Frames_results_correspondence = im2col(repmat(Frames_to_Seg, [Sub_Sample_Num,1]),[1 1]);
    Frames_results_correspondence = Frames_results_correspondence(1:nFrame);
    
    tic
    if(imageflattening_mode==2)
        %under mode 2, normalization is done without outlier-elimination
        % under mode 1, outliers are eliminated first
        
        % in mode 2( a collection of different images)
        
        % load one image to check the bit depth
        currentImg = movieData.channels_(iChannel).loadImage(1);
       
        % if input is double, use min 0(0.0000000001 for the log), 
        % max of uint8 if max is bigger then 1; else just set as 1
        
        if(isa(currentImg,'double'))
            if( flatten_method_ind==1)
                low_005_percentile = 0.0000000001;
            else
                low_005_percentile = 0;
            end
            
            if(max(max(currentImg))>1)
                high_995_percentile= 2^8-1;
            else
                high_995_percentile= 1;
            end
        end
          
        % set the max as the possible max, min as possible min 
        if( isa(currentImg,'uint16'))
            if( flatten_method_ind==1)
                low_005_percentile = 1;
            else
                low_005_percentile = 0;
            end
            % the reason for using 12 is 16bit images are usually dim if
            % opened in image viewer, since not 16 bits are fully used.
            % Checked again, should be a microscope issue, the image is 12
            % bits.
            if(max(max(currentImg))>2^12-1)
                high_995_percentile= 2^16-1;
            else                
                high_995_percentile= 2^12-1;                
            end
        end
        
        if( isa(currentImg,'uint8'))
            if( flatten_method_ind==1)
                low_005_percentile = 1;
            else
                low_005_percentile = 0;
            end
            high_995_percentile= 2^8-1;            
        end
        
        center_value=ones(movieData.nFrames_,1);        
        center_value_int=0;
        img_min=low_005_percentile;
        img_max=high_995_percentile;
        
    else
        % 0.5 percdentile and 99.5 percentile are used to eliminate the
        % outliers(extremely bright or dim spots)
        img_pixel_pool = [];
        for iFrame_subsample = 1 : length(Frames_to_Seg)
            iFrame = Frames_to_Seg(iFrame_subsample);
            currentImg = movieData.channels_(iChannel).loadImage(iFrame);
            
            if(size(currentImg,1)>100 || length(Frames_to_Seg)>20)
                smaller_currentImg = imresize(currentImg,[100 NaN]);
                if length(Frames_to_Seg)>50
                    smaller_currentImg = imresize(currentImg,[50 NaN]);
                end
            else
                smaller_currentImg = currentImg;
            end
            
            img_pixel_pool = [img_pixel_pool smaller_currentImg(:)];
        end
        [hist_all_frame, hist_bin] = hist(double(img_pixel_pool),100);
        
        img_pixel_pool = double(img_pixel_pool(:));
        nonzero_img_pixel_pool= img_pixel_pool(img_pixel_pool>0);
        
        low_005_percentile = prctile(img_pixel_pool,0.01);
        
        % for log mode, need to find a min bigger than 0
        if(low_005_percentile==0 && flatten_method_ind==1)
            low_005_percentile = prctile(img_pixel_pool,1);            
            if(low_005_percentile==0)
                low_005_percentile = 1*min(nonzero_img_pixel_pool);
            end
        end        
        
        % if not found the loop use 1 max
        high_995_percentile = prctile(img_pixel_pool,99.99);
        
        if exist('img_pixel_pool','var')
            clearvars img_pixel_pool;
        end
        
        img_min=low_005_percentile;
        img_max=high_995_percentile;
        
        currentImg_cell = cell(1,1);
        
        for iFrame_subsample = 1 : length(Frames_to_Seg)
            hist_this_frame = hist_all_frame(:,iFrame_subsample);
            ind = find(hist_this_frame==max(hist_this_frame));
            center_value(iFrame_subsample) = hist_bin(ind(1));
            if(ind(1)>1)
                center_value_m1(iFrame_subsample) = hist_bin(ind(1)-1);
            else
                center_value_m1(iFrame_subsample)= center_value(iFrame_subsample);
            end
            %       center_value(iFrame_subsample) = 1;
        end
        center_value_int = mean((center_value_m1+center_value)/2);
        
        % record the stat numbers
        % funParams.stat.low_005_percentile = low_005_percentile;
        % funParams.stat.high_995_percentile = high_995_percentile;
        % funParams.stat.center_value_int = center_value_int;
        % QZ I do not think those stat are saved somewhere -- start
        p.stat.low_005_percentile = low_005_percentile;
        p.stat.high_995_percentile = high_995_percentile;
        p.stat.center_value_int = center_value_int;
        % QZ -- end

        center_value = center_value/max(center_value);
        center_value = sqrt(center_value);
        center_value = imfilter(center_value,[1 2 3 9 3 2 1]/21,'replicate','same');
        center_value = center_value/max(center_value);
        
        %     center_value(:)=1;
%         img_min=0;
%         img_max=high_995_percentile- center_value_int;
%         
        img_min=low_005_percentile;
        img_max=high_995_percentile;
       
    end
    
    display(['Time for statistics of image intensity in Channel ',num2str(iChannel)]);
   
    toc
    
    %%
        
    % % Make output directory for the flattened images
    % ImageFlattenChannelOutputDir = movieData.processes_{indexFlattenProcess}.outFilePaths_{iChannel};
    % if (~exist(ImageFlattenChannelOutputDir,'dir'))
    %     mkdir(ImageFlattenChannelOutputDir);
    % end
    
    %%
    display('======================================');
   
    for iFrame_subsample = 1 : length(Frames_to_Seg)
        iFrame = Frames_to_Seg(iFrame_subsample);
        disp(['Frame: ',num2str(iFrame)]);
        
        tic
        
        % Read in the intensity image.
        currentImg = movieData.channels_(iChannel).loadImage(iFrame);
        currentImg = double(currentImg);
        
%         currentImg = currentImg - center_value_int;
        
        % Get rid of extreme noises
%         currentImg(find(currentImg>high_995_percentile- center_value_int))=high_995_percentile- center_value_int;
        currentImg(find(currentImg>high_995_percentile))= high_995_percentile;
        currentImg(find(currentImg<low_005_percentile))= low_005_percentile;
        
        % based on the given method index, do log or sqrt to flatten the image
        if flatten_method_ind == 1
            
            currentImg = log(currentImg);
            currentImg = currentImg - log(img_min);
            currentImg = currentImg/...
                (log((center_value(iFrame_subsample))*img_max) ...
                -log((center_value(iFrame_subsample))*img_min));
            
        else
            
            currentImg = currentImg - img_min;
            currentImg = currentImg/(center_value(iFrame_subsample))/(img_max- img_min);
            if flatten_method_ind == 2
                currentImg = (currentImg).^(1/2);
            end
            
            if flatten_method_ind == 3
                currentImg = (currentImg).^(2/3);
            end
        end
        
        
        % Smooth the image in requested
        if Gaussian_sigma > 0
            currentImg = imfilter(currentImg, fspecial('gaussian',round(5*Gaussian_sigma), Gaussian_sigma),'replicate','same');
        end
        
        currentImg_cell{iFrame}=currentImg;
        
        %% %tif stack cost too much memory, geta cell instead
        if(TimeFilterSigma > 0)
            for_temporal_filtering{iFrame_subsample} = currentImg;
        else
            % if no need for time filtering, save to disk
            for sub_i = 1 : Sub_Sample_Num
                if iFrame + sub_i-1 <= nFrame
                    
                    currentImg(currentImg<0)=0;
                    currentImg(currentImg>1)=1;
                    currentImg = currentImg*(2^8-1);
                    currentImg = uint8(currentImg);
                    
                    imwrite(currentImg,[ImageFlattenChannelOutputDir,filesep,'flatten_', ...
                        filename_short_strs{iFrame + sub_i-1},'.tif']);
                end
            end
        end
        
        toc
        
    end
    
    %% if temporal filtering is needed
    
    if(TimeFilterSigma > 0 && imageflattening_mode ==1)
        disp('Image Flattening in temporal filtering:');
        
        tic
        
        % Initialize the after filtering cell
        after_temporal_filtering = cell(1);
        
        % prepare the temporal filter
        FilterHalfLength = 2*ceil(TimeFilterSigma);
        temperal_filter = zeros(1,1,2*FilterHalfLength+1);
        H = fspecial('gaussian',2*FilterHalfLength+1, TimeFilterSigma);
        H_1D = H(FilterHalfLength+1,:);
        H_1D = H_1D/(sum(H_1D));
        temperal_filter(1,1,:) = H_1D(:);
        
        % cut the sequence into trunks(length pace), to reduce the need of stacking a
        % long long sequence into one 3D matrix which kills the memory
        
        pace = min(length(Frames_to_Seg),10);
        
        for iTfilter = 1 : pace : length(Frames_to_Seg)
            
            % these are the index of starting and ending of image to be
            % filtered, so could include the neighboring image on both ends
            start_iF_sub = max(1, iTfilter-ceil(TimeFilterSigma));
            end_iF_sub = min(length(Frames_to_Seg), iTfilter+pace-1+ceil(TimeFilterSigma));
            
            % this is the results, without the padding
            start_filtered = iTfilter;
            end_filtered = min(length(Frames_to_Seg), iTfilter+pace-1);
            
            % initialize the 3D tensor
            image_tensor_10 = zeros(size(currentImg,1),size(currentImg,2),end_iF_sub-start_iF_sub+1);
            
            for iF_sub = start_iF_sub : end_iF_sub
                image_tensor_10(:,:,  iF_sub - start_iF_sub +1) ...
                    = for_temporal_filtering{iF_sub};
            end
            
            % filter the 3D tensor
            time_filtered_10 = imfilter(image_tensor_10,temperal_filter,'replicate','same');
            
            % read out into the results, without padding
            for iF_sub = start_filtered : end_filtered
                after_temporal_filtering{iF_sub}=...
                    squeeze(time_filtered_10(:,:,iF_sub-start_iF_sub+1));
            end
            
        end
        
        % save the filtered to hard disk
        for iFrame_subsample = 1 : length(Frames_to_Seg)
            iFrame = Frames_to_Seg(iFrame_subsample);
            disp(['Image Flattening in temporal filtering, Frame: ',num2str(iFrame)]);
            
            currentImg = after_temporal_filtering{iFrame_subsample};
            currentImg(currentImg<0)=0;
            currentImg(currentImg>1)=1;
            
            currentImg = currentImg*(2^8-1);
            currentImg = uint8(currentImg);
            
            for sub_i = 1 : Sub_Sample_Num
                if iFrame + sub_i-1 <= nFrame
                    disp(['Frame: ',num2str(iFrame + sub_i-1)]);
                    
                    imwrite(currentImg,[ImageFlattenChannelOutputDir,filesep,'flatten_', ...
                        filename_short_strs{iFrame + sub_i-1},'.tif']);
                end
            end
        end
        disp('Total time in Image Flattening in filtering:');
        
        toc
    end
    
    %% if background removal is needed
    if background_removal_flag==1 % QZ this part is never called
        disp('Image Flattening in background removal:');
        
        tic
        % Background substraction for uneven illumination
        for iFrame_subsample = 1 : length(Frames_to_Seg)
            iFrame = Frames_to_Seg(iFrame_subsample);
            disp(['Image Flattening in back ground removal, Frame: ',num2str(iFrame)]);
            currentImg =  imread([ImageFlattenChannelOutputDir,filesep,'flatten_', ...
                filename_short_strs{iFrame + sub_i-1},'.tif']);
            
            I = double(currentImg);
            % Get the x and y of the surface
            [XI, YI] = meshgrid(1:size(I,2), 1:size(I,1));
            % Fit a polynomial to the surface x-y both up to 2nd order
            fit_sur = fit([YI(:),XI(:)],I(:), 'poly22', 'Robust', 'on');
            % Reconstract the surface for uneven background---without the first part (the DC part)
            Z_fit = fit_sur.p01.*XI+ fit_sur.p10.*YI +fit_sur.p11.*XI.*YI+ ...
                fit_sur.p20.*YI.*YI+fit_sur.p02.*XI.*XI;
            % Substract the background
            currentImg = I-Z_fit;
            currentImg = currentImg/255;% back to 0~1 for saving image tif
            
            currentImg(currentImg<0)=0;
            currentImg(currentImg>1)=1;
            
            currentImg = currentImg*(2^8-1);
            currentImg = uint8(currentImg);
            
            % Save to disk
            for sub_i = 1 : Sub_Sample_Num
                if iFrame + sub_i-1 <= nFrame
                    disp(['Frame: ',num2str(iFrame + sub_i-1)]);
                    
                    imwrite(currentImg,[ImageFlattenChannelOutputDir,filesep,'flatten_', ...
                        filename_short_strs{iFrame + sub_i-1},'.tif']);
                end
            end
            
        end
        
        disp('Total time in Image Flattening in background removal:');
        
        toc
    end
    
    %%
    % this the end of "for" of each channel
end
%%%% end of algorithm

fprintf('\n Finished Image Flatten Process! \n')

end