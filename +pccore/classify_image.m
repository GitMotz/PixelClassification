function [currentSegmentation,currentPSignal,currentPBackground] = classify_image(imName,hData,im,imSec,imTer,CellSeg)


%find the image case
imagecase = 'one_image';
if nargin == 3 % || ( nargin == 5 && isempty(imSec) && isempty(imTer) ) 
    imagecase = 'one_image';
elseif nargin == 4 % || ( nargin == 5 && ~isempty(imSec) && isempty(imTer) ) 
    imagecase = 'two_images';
elseif nargin == 5 %|| ( nargin == 5 && ~isempty(imSec) && ~isempty(imTer) ) 
    imagecase = 'three_images';
elseif nargin == 6
    if isempty(imSec) && isempty(imTer)
        imagecase = 'one_image';
    elseif ~isempty(imSec) && isempty(imTer)
        imagecase = 'two_images';
    elseif ~isempty(imSec) && ~isempty(imTer)
        imagecase = 'three_images';
    else
        imagecase = 'unknown';
    end    
else    
    imagecase = 'unknown';
end
 
%find the segmentation case
cellsegcase = 'cell_seg_false';
if nargin == 6 && ~isempty(CellSeg)
    cellsegcase = 'cell_seg_true';
end


% segment image
PadSize = 100;

fprintf('%s: getting features of first image.\n',mfilename);
tim = padarray(im,[PadSize PadSize],'symmetric','both');
pixelFeatures = pcsf.GeneratePixelFeatures(tim,hData.structSettings.FeatureList);
pixelFeatures = pcsf.getFeatureSection2(pixelFeatures,PadSize);

if isfield(hData.structSettings,'UsePreSegmentation')
    
    [pathstr,namestr] = fileparts(imName);
    pathstr = fullfile(pathstr(1:end-5),'SEGMENTATION');
    namestr = [namestr(1:end-9) 'A01Z01C01_Segmented' hData.structSettings.UsePreSegmentation '.png'];
    seg_image = imread(fullfile(pathstr,namestr));
    [d1,d2,d3,norm_int] = pcsf.compute_distance_features(im,seg_image);
    pixelFeatures{end+1}=d1;
    pixelFeatures{end+1}=d2;
    pixelFeatures{end+1}=d3;
    pixelFeatures{end+1}=norm_int;

end



if  strcmp(imagecase,'two_images') ||  strcmp(imagecase,'three_images')
    fprintf('%s: getting features of second image.\n',mfilename);
    tim = padarray(imSec,[PadSize PadSize],'symmetric','both');
    tf = pcsf.GeneratePixelFeatures(tim,hData.structSettings.SecondaryFeatureList);
    tf = pcsf.getFeatureSection2(tf,PadSize);
    for iFeat = 1:length(tf)
        pixelFeatures{end+1}=tf{iFeat};
    end
    clear tf
end

if  strcmp(imagecase,'three_images')
    fprintf('%s: getting features of third image.\n',mfilename);
    tim = padarray(imTer,[PadSize PadSize],'symmetric','both');
    tf = pcsf.GeneratePixelFeatures(tim,hData.structSettings.TertiaryFeatureList);
    tf = pcsf.getFeatureSection2(tf,PadSize);
    for iFeat = 1:length(tf)
        pixelFeatures{end+1}=tf{iFeat};
    end
    clear tf
end


if isfield(hData.structSettings,'UseProbabilityMask')    
    [pathstr,namestr] = fileparts(imName);
    pathstr = fullfile(pathstr(1:end-5),'PreSegTest');
    namestr = [namestr '_Probability' hData.structSettings.UseProbabilityMask '.png'];
    prob_image = double(imread(fullfile(pathstr,namestr)))./(2^16);
    pixelFeatures{end+1}=prob_image;
end



if isfield(hData.structSettings,'UseSegmentationCellFeatures')
    [pathstr,namestr] = fileparts(imName);
    pathstr = fullfile(pathstr(1:end-5),'SEGMENTATION');
    namestr = [namestr(1:end-9) 'A01Z01C01_Segmented' hData.structSettings.UseSegmentationCellFeatures '.png'];
    seg_image = imread(fullfile(pathstr,namestr));
    
    ix1 = unique(seg_image(:));
    ix1(ix1==0)=[];
    
    m1 = zeros(size(im));
    m2 = zeros(size(im));


    for iOb = 1:length(ix1)
        f = seg_image==ix1(iOb);
        d = im(f(:));
        m1(f(:)) = kurtosis(d);
        m2(f(:)) = skewness(d);
       
    end

    pixelFeatures{end+1}=m1;
    pixelFeatures{end+1}=m2;
end






if  strcmp(imagecase,'unknown')
    error('The image case is unknown.')
end


switch cellsegcase 
    case 'cell_seg_true'
        
        if  strcmp(imagecase,'two_images') ||  strcmp(imagecase,'three_images')
            boundaries = round(linspace(1,size(pixelFeatures{1},2)+1,7));
            for iIm = 1:length(boundaries)-1
               ix1 = boundaries(iIm);
               ix2 = boundaries(iIm+1)-1;
               for iFeat = 1:length(pixelFeatures)
                    tempPixelFeatures{iFeat} = pixelFeatures{iFeat}(:,ix1:ix2);
               end
               [tempCurrentPSignal{iIm},tempCurrentPBackground{iIm}] = pccore.basic_classification(hData,tempPixelFeatures,CellSeg(:,ix1:ix2));
            end
            currentPSignal = cat(2,tempCurrentPSignal{:});
            currentPBackground = cat(2,tempCurrentPBackground{:});
        else
            [currentPSignal,currentPBackground] = pccore.basic_classification(hData,pixelFeatures,CellSeg);
        end
        
        
    case 'cell_seg_false';
        
        if  strcmp(imagecase,'two_images') ||  strcmp(imagecase,'three_images')
            boundaries = round(linspace(1,size(pixelFeatures{1},2)+1,7));
            for iIm = 1:length(boundaries)-1
               ix1 = boundaries(iIm);
               ix2 = boundaries(iIm+1)-1;
               for iFeat = 1:length(pixelFeatures)
                    tempPixelFeatures{iFeat} = pixelFeatures{iFeat}(:,ix1:ix2);
               end
               [tempCurrentPSignal{iIm},tempCurrentPBackground{iIm}] = pccore.basic_classification(hData,tempPixelFeatures);
            end
            currentPSignal = cat(2,tempCurrentPSignal{:});
            currentPBackground = cat(2,tempCurrentPBackground{:});
        else
            [currentPSignal,currentPBackground] = pccore.basic_classification(hData,pixelFeatures);
        end

        
        
        
end




[currentPSignal,currentPBackground] = pccore.smooth_scores(hData,currentPSignal,currentPBackground);

%currentPSignal = Score_1;
%currentPBackground = Score_2;
currentSegmentation = currentPSignal>hData.svm.threshold;
currentSegmentation = imfill(bwlabel(currentSegmentation));

