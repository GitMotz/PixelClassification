% this script runs the classification model on all images of the example
% data set

pathToDir = pwd;
% pathToModel = fullfile(pathToDir,'\ExampleDataSet\PCProject01.mat'); % (change name of model)
pathToModel = fullfile(pathToDir,'\ExampleDataSet\segmented_splicing_speckles\PixelClassification_SplicingSpeckles.mat');

% create settings file for function run_pix_class
StructSettings.classification_file = pathToModel;
StructSettings.RootPathBrutus = fullfile(pathToDir,'\ExampleDataSet\');
StructSettings.RootPath = fullfile(pathToDir,'\ExampleDataSet\');
StructSettings.ObjectName = 'splicing_speckles';
StructSettings.IllumCorrFile = 'Measurements_batch_illcor_channel003_zstack000.mat';
ChannelId = {'C03.png'};            
strDir = CPdir(fullfile(StructSettings.RootPathBrutus,'TIFF'));
strDir = {strDir(:).name};
f = ~cellfun(@isempty,strfind(strDir,ChannelId{1}));
cellImages =  strDir(f);
StructSettings.jobs.ImageList = cellImages;

save(fullfile(pathToDir,'\ExampleDataSet\run_pix_class_settings.mat'),'StructSettings','-v7.3' )


% apply pixel classification to all images
run_pix_class(fullfile(pathToDir,'\ExampleDataSet\run_pix_class_settings.mat'));

