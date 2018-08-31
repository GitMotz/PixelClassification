function run_pix_class(JobSettingsFile)

% load the job settings file
StructSettings = importdata(JobSettingsFile);

imageList = StructSettings.jobs.ImageList;
classModel = importdata(StructSettings.classification_file);

% load illumination correction files
corrName = fullfile(StructSettings.RootPathBrutus,'TIFF',StructSettings.IllumCorrFile);
stat_values = importdata(corrName);

parfor i = 1:size(imageList,2)
    
    % load image
    imName = fullfile(StructSettings.RootPathBrutus,'TIFF',StructSettings.jobs.ImageList{i});
    im = double(imread(imName));
    im = IllumCorrect(im,stat_values.mean,stat_values.std,1);
    
    
    % get nuclear segmentation
    try
        fprintf('%s: Trying to get the nuclear segmentation.\n',mfilename)
        cImName = StructSettings.jobs.ImageList{i};
        expression = 'A(\d+)Z';
        cImName = regexprep(cImName,expression,'A01Z','preservecase');
        expression = 'C(\d+).png';
        cImName = regexprep(cImName,expression,'C01_SegmentedNuclei.png');
        cImName = fullfile(StructSettings.RootPathBrutus,'SEGMENTATION',cImName);
        CellSeg = imread(cImName);
        cellsegcase = 'cell_seg_true';
    catch
        warning('No nuclear segmentation found.')
        cellsegcase = 'cell_seg_false';
    end
    
    % apply model 
    switch cellsegcase
        case 'cell_seg_false'
            [currentSegmentation,~] = pccore.classify_image(imName,classModel,im);
        case 'cell_seg_true'
            [currentSegmentation,~] = pccore.classify_image(imName,classModel,im,[],[],CellSeg);
    end
    
   
    if any(currentSegmentation(:)>0)
        IM1 = currentSegmentation;
        [CurrentObjNhood,CurrentObjLabels] = bwdist(IM1 );
        CurrentObjNhood = double(CurrentObjNhood < 3) .* double(IM1(CurrentObjLabels));
        
        IM2 = imfill(CurrentObjNhood);
        
        BW2 = zeros(size(IM2));
        uniqueSth  = unique(IM2(:));
        uniqueSth(uniqueSth==0)=[];
        for iObject = 1:length(uniqueSth)
            ix = bwmorph(IM2==uniqueSth(iObject),'remove');
            BW2(ix(:)) = 1;
        end
        
        IM3=IM2;
        IM3(BW2>0)=0;
        IM3 = imfill(IM3);
        IM3 = bwlabel(IM3);
        
        se = strel('disk',4);
        IM3 = (imerode(IM3,se));
        
        IM4 = (IM3+IM1)>0;
        IM4 = bwlabel(IM4);
        
        currentSegmentation = IM4;
           
    end
    
    % write segmentation image
    OutputName = [StructSettings.jobs.ImageList{i}(1:end-4) '_Segmented_' StructSettings.ObjectName '.png'];
    strFullFileName = fullfile(StructSettings.RootPathBrutus,'TEST_seg',OutputName);
    imwrite(uint16(currentSegmentation),strFullFileName,'png','BitDepth',16);
    fprintf('Wrote image %s.\n',StructSettings.jobs.ImageList{i})
    
    
end

end



