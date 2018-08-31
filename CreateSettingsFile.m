
% this script generates a setting file that is required for ClassifyPixels.m

pathToDir = pwd;
pathToImages = fullfile(pathToDir,'\ExampleDataSet\TIFF');
pathOutput = fullfile(pathToDir,'\settings.mat');

% create paths to images
ChannelId = {'C03.png'};
expression = 'C0(\d).png';
intChannelNumber = regexp(ChannelId,expression,'tokens');

structSettings=struct();

fileDir = dir(pathToImages);
allFiles = {fileDir(:).name};
ix_imageFiles = ~cellfun(@isempty,strfind(allFiles,ChannelId{1}));
imageFiles = allFiles(ix_imageFiles);

% write paths to images
structSettings.pathToImages = {};
for i = 1:size(imageFiles,2)
    structSettings.pathToImages{i} = fullfile(pathToImages,imageFiles{i});
end

% write paths to illumination correction file
structSettings.pathToIllumCorr = {};
for i = 1:size(imageFiles,2)
    structSettings.pathToIllumCorr{i} = fullfile(pathToImages,sprintf('Measurements_batch_illcor_channel00%s_zstack000.mat',intChannelNumber{1}{1}{1}));
end

% write pixel features to be calculated
structSettings.FeatureList = {'log_5_2';'log_8_4';'log_20_10';'gau_20_10';'gau_9_4';'gau_5_2';'diffgau_1';'diffgau_2';...
    'average_5';'average_10';'average_20';'edges';'hole_log_8_4';'morph_opening';'morph_close';'deconv_10_7';...
    'deconv_10_3';'std_dif_1';'zscore_1'};

save(pathOutput,'structSettings','-v7.3');


% All available pixel features (see also GeneratePixelFeatures.m):
% {'log_5_2';'log_5_3';'log_6_3';'log_8_4';'log_10_5';'log_20_10';
%  'gau_20_10';'gau_9_4';'gau_5_2';'diffgau_1';'diffgau_2';
%  'average_3';'average_5';'average_10';'average_20';
%  'edges';'simple_edges';'edge_log_8_4';'weird_log_8_4';'entropy_log_8_4';'entropy_log_5_2';
%  'morph_opening';'morph_close';'morph_close_small';'loc_back_small';'loc_back_large';
%  'deconv_10_7';'deconv_10_3';'deconv_5_3';'std_dif_1';'zscore_1';
%  'watershed_lines';'watershed_lines_smooth';'watershed_lines_smooth_low';'watershed_lines_smooth_very_low';...
%  'tophat_12';'orig_tophat_30'}; 
