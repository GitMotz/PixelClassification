function varargout = ClassifyPixels(varargin)
% CLASSIFYPIXELS MATLAB code for ClassifyPixels.fig
%      CLASSIFYPIXELS, by itself, creates a new CLASSIFYPIXELS or raises the existing
%      singleton*.
%
%      H = CLASSIFYPIXELS returns the handle to a new CLASSIFYPIXELS or the handle to
%      the existing singleton*.
%
%      CLASSIFYPIXELS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLASSIFYPIXELS.M with the given input arguments.
%
%      CLASSIFYPIXELS('Property','Value',...) creates a new CLASSIFYPIXELS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ClassifyPixels_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ClassifyPixels_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ClassifyPixels

% Last Modified by GUIDE v2.5 25-Apr-2016 17:49:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ClassifyPixels_OpeningFcn, ...
    'gui_OutputFcn',  @ClassifyPixels_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



% --- Executes just before ClassifyPixels is made visible.
function ClassifyPixels_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ClassifyPixels (see VARARGIN)

% Choose default command line output for ClassifyPixels
handles.output = hObject;

% Load icon
iconPixel = imread('iCON.png');
iconPixel(iconPixel>=150)=max(iconPixel(:));
iconPixel(iconPixel<150)=0;
iconPixel = iconPixel(31:end-30,31:end-30,:);
imshow(iconPixel, 'Parent', handles.axes3);

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes ClassifyPixels wait for user response (see UIRESUME)
% uiwait(handles.classify_pixels_master);


% --- Outputs from this function are returned to the command line.
function varargout = ClassifyPixels_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in class_1.
function class_1_Callback(hObject, eventdata, handles)
% hObject    handle to class_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hData
hData.svm.currentMethod = 'Signal';
set(handles.current_color,...
    'BackgroundColor',[0 1 0],...
    'String','Signal');

% --- Executes on button press in class_2.
function class_2_Callback(hObject, eventdata, handles)
% hObject    handle to class_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hData
hData.svm.currentMethod = 'Background';
set(handles.current_color,...
    'BackgroundColor',[1 0 0],...
    'String','Background');

% --- Executes on button press in eraser.
function eraser_Callback(hObject, eventdata, handles)
% hObject    handle to eraser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hData
hData.svm.currentMethod = 'Eraser';
set(handles.current_color,...
    'BackgroundColor',[0 0 0],...
    'String','Eraser');

% --- Executes on button press in next_image.
function next_image_Callback(hObject, eventdata, handles)
% hObject    handle to next_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hData
[hData.svm.XTraining,ia]=unique([hData.svm.XTraining;hData.svm.currentXTraining],'rows');
hData.svm.YTraining=[hData.svm.YTraining;hData.svm.currentYTraining];
hData.svm.YTraining=hData.svm.YTraining(ia);

hData.currentSectionIndex = 1;
hData.currentImageIndex = hData.currentImageIndex+1;
fprintf('%s: The current image index is: %.4d\n',mfilename,hData.currentImageIndex)
if hData.currentImageIndex>length(hData.structSettings.pathToImages)
    hData.currentImageIndex = 1;
end


randSection = randperm(25,1);
hData.currentSectionIndex = randSection;
set(handles.Section_Number,'String', num2str(hData.currentSectionIndex));

PadSize = 100;

disp('Loading next image and computing features, please wait...')
hData.currentOrigImage = double(imread(hData.structSettings.pathToImages{hData.currentImageIndex}));
load(hData.structSettings.pathToIllumCorr{hData.currentImageIndex});
hData.currentOrigImageUncroped = IllumCorrect(hData.currentOrigImage,stat_values.mean,stat_values.std,1);
hData.currentOrigImageBkup = hData.currentOrigImage;
hData.currentOrigImage = pcsf.getPreImageSection(hData.currentOrigImageUncroped,PadSize,hData.currentSectionIndex);
hData.pixelFeatures = pcsf.GeneratePixelFeatures(hData.currentOrigImage,hData.structSettings.FeatureList);

if isfield(hData.structSettings,'UsePreSegmentation')
    [pathstr,namestr] = fileparts(hData.structSettings.pathToImages{hData.currentImageIndex});
    pathstr = fullfile(pathstr(1:end-5),'SEGMENTATION');
    namestr = [namestr(1:end-9) 'A01Z01C01_Segmented' hData.structSettings.UsePreSegmentation '.png'];
    seg_image = imread(fullfile(pathstr,namestr));
    seg_image = pcsf.getPreImageSection(seg_image,PadSize,hData.currentSectionIndex);
    hData.currentSegImage = seg_image;
    hData.currentSegSection = pcsf.getImageSection2(hData.currentSegImage,PadSize);

    [d1,d2,d3,norm_int] = pcsf.compute_distance_features(hData.currentOrigImage,seg_image);
    hData.pixelFeatures{end+1}=d1;
    hData.pixelFeatures{end+1}=d2;
    hData.pixelFeatures{end+1}=d3;
    hData.pixelFeatures{end+1}=norm_int;
    
else
    [pathstr,namestr] = fileparts(hData.structSettings.pathToImages{hData.currentImageIndex});
    pathstr = fullfile(pathstr(1:end-5),'SEGMENTATION');
    namestr = [namestr(1:end-9) 'A01Z01C01_SegmentedNuclei.png'];
    seg_image = imread(fullfile(pathstr,namestr));
    seg_image = pcsf.getPreImageSection(seg_image,PadSize,hData.currentSectionIndex);
    hData.currentSegImage = seg_image;
    hData.currentSegSection = pcsf.getImageSection2(hData.currentSegImage,PadSize);
end
    
    


% load secondary image if possible
if isfield(hData.structSettings,'pathToSecondaryImages')
    hData.currentSecOrigImage = double(imread(hData.structSettings.pathToSecondaryImages{hData.currentImageIndex}));
    load(hData.structSettings.pathToSecondaryIllumCorr{hData.currentImageIndex});
    hData.currentSecOrigImageUncroped  = IllumCorrect(hData.currentSecOrigImage,stat_values.mean,stat_values.std,1);
    hData.currentSecOrigImage = pcsf.getPreImageSection(hData.currentSecOrigImageUncroped,PadSize,hData.currentSectionIndex);
    tf = pcsf.GeneratePixelFeatures(hData.currentSecOrigImage,hData.structSettings.SecondaryFeatureList);
    for iFeat = 1:length(tf)
        hData.pixelFeatures{end+1}=tf{iFeat};
    end
end

% load secondary image if possible
if isfield(hData.structSettings,'pathToTertiaryImages')
    hData.currentTerOrigImage = double(imread(hData.structSettings.pathToTertiaryImages{hData.currentImageIndex}));
    load(hData.structSettings.pathToTertiaryIllumCorr{hData.currentImageIndex});
    hData.currentTerOrigImageUncroped = IllumCorrect(hData.currentTerOrigImage,stat_values.mean,stat_values.std,1);
    hData.currentTerOrigImage = pcsf.getPreImageSection(hData.currentTerOrigImageUncroped,PadSize,hData.currentSectionIndex);
    tf = pcsf.GeneratePixelFeatures(hData.currentTerOrigImage,hData.structSettings.TertiaryFeatureList);
    for iFeat = 1:length(tf)
        hData.pixelFeatures{end+1}=tf{iFeat};
    end
end


if isfield(hData.structSettings,'UseProbabilityMask')
    [pathstr,namestr] = fileparts(hData.structSettings.pathToImages{hData.currentImageIndex});
    pathstr = fullfile(pathstr(1:end-5),'PreSegTest');
    namestr = [namestr '_Probability' hData.structSettings.UseProbabilityMask '.png'];
    prob_image = double(imread(fullfile(pathstr,namestr)))./(2^16);
    prob_image = pcsf.getPreImageSection(prob_image,PadSize,hData.currentSectionIndex);
    hData.pixelFeatures{end+1}=prob_image;
end



if isfield(hData.structSettings,'UseSegmentationCellFeatures')
    [pathstr,namestr] = fileparts(hData.structSettings.pathToImages{hData.currentImageIndex});
    pathstr = fullfile(pathstr(1:end-5),'SEGMENTATION');
    namestr = [namestr(1:end-9) 'A01Z01C01_Segmented' hData.structSettings.UseSegmentationCellFeatures '.png'];
    seg_image = imread(fullfile(pathstr,namestr));
    seg_image = pcsf.getPreImageSection(seg_image,PadSize,hData.currentSectionIndex);
    
    ix1 = unique(seg_image(:));
    ix1(ix1==0)=[];
    
    m1 = zeros(size(hData.currentOrigImage));
    m2 = zeros(size(hData.currentOrigImage));


    for iOb = 1:length(ix1)
        f = seg_image==ix1(iOb);
        d = hData.currentOrigImage(f(:));
        m1(f(:)) = kurtosis(d);
        m2(f(:)) = skewness(d);
       
    end

    hData.pixelFeatures{end+1}=m1;
    hData.pixelFeatures{end+1}=m2;
end






hData.currentOrigSection = pcsf.getImageSection2(hData.currentOrigImage,PadSize);
hData.currentSectionPixelFeatures = pcsf.getFeatureSection2(hData.pixelFeatures,PadSize);
disp('Loading done, thanks for waiting.')

hData.lb = str2double( get(handles.LowQuantile,'String') );
hData.ub = str2double( get(handles.HighQuantile,'String') );
lb = quantile(hData.currentOrigImage(:),hData.lb);
ub = quantile(hData.currentOrigImage(:),hData.ub);
%hData.ub = quantile(hData.currentOrigImage(:),0.99995);
hData.currentOrigVisualizaion = hData.currentOrigSection-lb;
hData.currentOrigVisualizaion(hData.currentOrigVisualizaion<0)=0;
hData.currentOrigVisualizaion(hData.currentOrigVisualizaion>ub)=ub;

%hData.currentOrigVisualizaion = rescale_image(hData.currentOrigSection,hData.lb,hData.ub);
hData.svm.currentBackgrLabels = false(size(hData.currentOrigVisualizaion));
hData.svm.currentSignalLabels = false(size(hData.currentOrigVisualizaion));

R=hData.currentOrigVisualizaion./max(hData.currentOrigVisualizaion(:));
R(edge(hData.currentSegSection))=0;
G=hData.currentOrigVisualizaion./max(hData.currentOrigVisualizaion(:));
G(edge(hData.currentSegSection))=1;
B=hData.currentOrigVisualizaion./max(hData.currentOrigVisualizaion(:));
B(edge(hData.currentSegSection))=1;


hData.currentRGB=[];
hData.currentRGB(:,:,1) = R;
hData.currentRGB(:,:,2) = G;
hData.currentRGB(:,:,3) = B;
hData.currentRGB=hData.currentRGB./max(hData.currentRGB(:));

set(handles.axes1);
imagesc(hData.currentRGB, 'Parent', handles.axes1);
colormap('gray')
disp('Image:')
disp(hData.structSettings.pathToImages{hData.currentImageIndex})

% --- Executes on button press in next_section.
function next_section_Callback(hObject, eventdata, handles)
% hObject    handle to next_section (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hData

hData.currentSectionIndex = str2double(get(handles.Section_Number,'String'));
if hData.currentSectionIndex>25
    warning('The maximum section number is 25.\n');
    hData.currentSectionIndex = 25;
    set(handles.Section_Number,'String', num2str(hData.currentSectionIndex))   
end

PadSize = 100;

disp('Loading next image and computing features, please wait...')
hData.currentOrigImage = pcsf.getPreImageSection(hData.currentOrigImageUncroped,PadSize,hData.currentSectionIndex);
hData.pixelFeatures = pcsf.GeneratePixelFeatures(hData.currentOrigImage,hData.structSettings.FeatureList);


if isfield(hData.structSettings,'UsePreSegmentation')
   [pathstr,namestr] = fileparts(hData.structSettings.pathToImages{hData.currentImageIndex});
    pathstr = fullfile(pathstr(1:end-5),'SEGMENTATION');
    namestr = [namestr(1:end-9) 'A01Z01C01_Segmented' hData.structSettings.UsePreSegmentation '.png'];
    seg_image = imread(fullfile(pathstr,namestr));
    seg_image = pcsf.getPreImageSection(seg_image,PadSize,hData.currentSectionIndex);
    hData.currentSegImage = seg_image;
    hData.currentSegSection = pcsf.getImageSection2(hData.currentSegImage,PadSize);
    
    [d1,d2,d3,norm_int] = pcsf.compute_distance_features(hData.currentOrigImage,seg_image);
    hData.pixelFeatures{end+1}=d1;
    hData.pixelFeatures{end+1}=d2;
    hData.pixelFeatures{end+1}=d3;
    hData.pixelFeatures{end+1}=norm_int;

else
    [pathstr,namestr] = fileparts(hData.structSettings.pathToImages{hData.currentImageIndex});
    pathstr = fullfile(pathstr(1:end-5),'SEGMENTATION');
    namestr = [namestr(1:end-9) 'A01Z01C01_SegmentedNuclei.png'];
    seg_image = imread(fullfile(pathstr,namestr));
    seg_image = pcsf.getPreImageSection(seg_image,PadSize,hData.currentSectionIndex);
    hData.currentSegImage = seg_image;
    hData.currentSegSection = pcsf.getImageSection2(hData.currentSegImage,PadSize);
end


% load secondary image if possible
if isfield(hData.structSettings,'pathToSecondaryImages')
    hData.currentSecOrigImage = pcsf.getPreImageSection(hData.currentSecOrigImageUncroped,PadSize,hData.currentSectionIndex);
    tf = pcsf.GeneratePixelFeatures(hData.currentSecOrigImage,hData.structSettings.SecondaryFeatureList);
    for iFeat = 1:length(tf)
        hData.pixelFeatures{end+1}=tf{iFeat};
    end
end

% load secondary image if possible
if isfield(hData.structSettings,'pathToTertiaryImages')
    hData.currentTerOrigImage = pcsf.getPreImageSection(hData.currentTerOrigImageUncroped,PadSize,hData.currentSectionIndex);
    tf = pcsf.GeneratePixelFeatures(hData.currentTerOrigImage,hData.structSettings.TertiaryFeatureList);
    for iFeat = 1:length(tf)
        hData.pixelFeatures{end+1}=tf{iFeat};
    end
end

if isfield(hData.structSettings,'UseProbabilityMask')
    [pathstr,namestr] = fileparts(hData.structSettings.pathToImages{hData.currentImageIndex});
    pathstr = fullfile(pathstr(1:end-5),'PreSegTest');
    namestr = [namestr '_Probability' hData.structSettings.UseProbabilityMask '.png'];
    prob_image = double(imread(fullfile(pathstr,namestr)))./(2^16);
    prob_image = pcsf.getPreImageSection(prob_image,PadSize,hData.currentSectionIndex);
    hData.pixelFeatures{end+1}=prob_image;
end

if isfield(hData.structSettings,'UseSegmentationCellFeatures')
    [pathstr,namestr] = fileparts(hData.structSettings.pathToImages{hData.currentImageIndex});
    pathstr = fullfile(pathstr(1:end-5),'SEGMENTATION');
    namestr = [namestr(1:end-9) 'A01Z01C01_Segmented' hData.structSettings.UseSegmentationCellFeatures '.png'];
    seg_image = imread(fullfile(pathstr,namestr));
    seg_image = pcsf.getPreImageSection(seg_image,PadSize,hData.currentSectionIndex);
    
    ix1 = unique(seg_image(:));
    ix1(ix1==0)=[];
    
    m1 = zeros(size(hData.currentOrigImage));
    m2 = zeros(size(hData.currentOrigImage));


    for iOb = 1:length(ix1)
        f = seg_image==ix1(iOb);
        d = hData.currentOrigImage(f(:));
        m1(f(:)) = kurtosis(d);
        m2(f(:)) = skewness(d);
       
    end

    hData.pixelFeatures{end+1}=m1;
    hData.pixelFeatures{end+1}=m2;
end



hData.currentOrigSection = pcsf.getImageSection2(hData.currentOrigImage,PadSize);
hData.currentSectionPixelFeatures = pcsf.getFeatureSection2(hData.pixelFeatures,PadSize);
disp('Loading done, thanks for waiting.')

hData.lb = str2double( get(handles.LowQuantile,'String') );
hData.ub = str2double( get(handles.HighQuantile,'String') );
lb = quantile(hData.currentOrigImage(:),hData.lb);
ub = quantile(hData.currentOrigImage(:),hData.ub);
%hData.ub = quantile(hData.currentOrigImage(:),0.99995);
hData.currentOrigVisualizaion = hData.currentOrigSection-lb;
hData.currentOrigVisualizaion(hData.currentOrigVisualizaion<0)=0;
hData.currentOrigVisualizaion(hData.currentOrigVisualizaion>ub)=ub;

%hData.currentOrigVisualizaion = rescale_image(hData.currentOrigSection,hData.lb,hData.ub);

hData.svm.currentBackgrLabels = false(size(hData.currentOrigVisualizaion));
hData.svm.currentSignalLabels = false(size(hData.currentOrigVisualizaion));


R=hData.currentOrigVisualizaion./max(hData.currentOrigVisualizaion(:));
R(edge(hData.currentSegSection))=0;
G=hData.currentOrigVisualizaion./max(hData.currentOrigVisualizaion(:));
G(edge(hData.currentSegSection))=1;
B=hData.currentOrigVisualizaion./max(hData.currentOrigVisualizaion(:));
B(edge(hData.currentSegSection))=1;


hData.currentRGB=[];
hData.currentRGB(:,:,1) = R;
hData.currentRGB(:,:,2) = G;
hData.currentRGB(:,:,3) = B;
hData.currentRGB=hData.currentRGB./max(hData.currentRGB(:));

set(handles.axes1);
imagesc(hData.currentRGB, 'Parent', handles.axes1);
colormap('gray')
disp('Image:')
disp(hData.structSettings.pathToImages{hData.currentImageIndex})










% --- Executes on button press in train_svm.
function train_svm_Callback(hObject, eventdata, handles)
% hObject    handle to train_svm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hData



svm_type_val = get(handles.svm_type,'Value');
svm_type_text = get(handles.svm_type,'String');
hData.svm.kernel = svm_type_text{svm_type_val};
hData.svm.svm_type_val = svm_type_val;
[hData.svm.XTraining,ia]=unique([hData.svm.XTraining;hData.svm.currentXTraining],'rows');
hData.svm.YTraining=[hData.svm.YTraining;hData.svm.currentYTraining];
hData.svm.YTraining=hData.svm.YTraining(ia);

XTraining = hData.svm.XTraining;

% make sure they are single values
hData.svm.XTraining = single(hData.svm.XTraining);


YTraining = hData.svm.YTraining;
kernel = hData.svm.kernel; 
if isfield(hData.svm,'best_weights');
    for i = 1:size(XTraining,2)
        XTraining(:,i)=XTraining(:,i).*hData.svm.best_weights(i);
    end
end


[hData.svm.SVMStruct,hData.svm.CompactSVMModel,hData.svm.oberror] = pccore.svm_train(...
    XTraining,YTraining,kernel);

% --- Executes on button press in classify_svm.
function classify_svm_Callback(hObject, eventdata, handles)
% hObject    handle to classify_svm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hData
segmentation=zeros(size(hData.currentOrigSection));
Score_1=segmentation;
Score_2=segmentation;
segmentation(:)=1:length(segmentation(:));

% build the data for classification of the curret segment
allData = [];
for i = 1:length(hData.currentSectionPixelFeatures)
    allData = [allData hData.currentSectionPixelFeatures{i}(:)];
end

if isfield(hData.svm,'best_weights');
    for i = 1:size(allData,2)
        allData(:,i)=allData(:,i).*hData.svm.best_weights(i);
    end
end

allData = single(allData);

% for i = 1:size(allData,2)
%      allData(:,i)=(allData(:,i)-nanmean(hData.svm.XTraining(:,i)))./nanstd(hData.svm.XTraining(:,i));
% end
% classify image by collumn to reduce memory requirement
fprintf('%s: Classifying current section. This will take a moment.\n',mfilename);
tic
n = 1;
switch hData.svm.kernel
    case 'GentleBoost'
        [a, b]  = predict(hData.svm.CompactSVMModel, allData);
        Score_1(:)  = b(:,1);
        Score_2(:)  = b(:,2);
        
    case 'SVMEnsemble'
        
        for i = segmentation
            %[~, p]  = predict(hData.svm.CompactSVMModel1, allData(i,:));
            [a, b]  = pccore.svm_ensemble_predict(hData.svm.CompactSVMModel, [allData(i,:)]);
            segmentation(i) = a(:,:,1);
            b = nanmean(b,3);
            Score_1(i)  = b(:,1);
            Score_2(i)  = b(:,2);
            n = 1+n
        end
        
        
    otherwise
        for i = segmentation
            %[~, p]  = predict(hData.svm.CompactSVMModel1, allData(i,:));
            [a, b]  = predict(hData.svm.CompactSVMModel, [allData(i,:)]);
            try
                segmentation(i) = a;
            catch
                segmentation(i) = cellfun(@str2double,a);
            end
            Score_1(i)  = b(:,1);
            Score_2(i)  = b(:,2);
            n = 1+n;
        end
        
end

clear allData
fprintf('%s: Classification Finished. ',mfilename);
toc

% perform and save sermentation
hData.svm.smooth_p = get(handles.SmoothProb,'Value');
hData.svm.smooth_pm = get(handles.SmoothMore,'Value');
hData.svm.smooth_p_param = []; 
if hData.svm.smooth_p
    if hData.svm.smooth_pm
        h = fspecial('gaussian', 5, 1.5); h = h./sum(h(:));
        Score_1 = imfilter(Score_1,h,'symmetric');
        Score_2 = imfilter(Score_2,h,'symmetric');
        hData.svm.smooth_p_param = [5 1.5];
    else
        h = fspecial('gaussian', 5, 0.75); h = h./sum(h(:));
        Score_1 = imfilter(Score_1,h,'symmetric');
        Score_2 = imfilter(Score_2,h,'symmetric');    
        hData.svm.smooth_p_param = [5 0.7];
    end
else
    h = fspecial('gaussian', 3, 0.5); h = h./sum(h(:));
    Score_1 = imfilter(Score_1,h,'symmetric');
    Score_2 = imfilter(Score_2,h,'symmetric');
    hData.svm.smooth_p_param = [3 0.5];
end    
    
switch hData.svm.kernel
    case 'GentleBoost'
    otherwise
        Score_1=Score_1./(Score_1+Score_2);
        Score_2=Score_2./(Score_1+Score_2);
end

hData.svm.currentPSignal = Score_1;
hData.svm.currentPBackground = Score_2;
hData.svm.threshold = str2double(get(handles.seg_threshold,'string'));
hData.svm.currentSegmentation = hData.svm.currentPSignal>hData.svm.threshold;
hData.svm.currentSegmentation = imfill(bwlabel(hData.svm.currentSegmentation));
hData.svm.currentSegmentation = hData.svm.currentSegmentation > 0;


% --- Executes on button press in show_original_image.
function show_original_image_Callback(hObject, eventdata, handles)
% hObject    handle to show_original_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hData
fprintf('%s: Displaying original image as a separate figure.\n',mfilename)
im = rescale_image(hData.currentOrigImageBkup,hData.lb,hData.ub);
im = im./max(im(:));
im(:,:,2)=im;
im(:,:,3)=im(:,:,1);

showOrigSeg = get(handles.show_orig_im_segmentation,'Value');
if showOrigSeg
    im = addCellNucSeg(im);
end
figure; imagesc(im)

% --- Executes on button press in show_section_image.
function show_section_image_Callback(hObject, eventdata, handles)
% hObject    handle to show_section_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hData
set(handles.axes1);
imagesc(hData.currentRGB,'Parent',handles.axes1);

% --- Executes on button press in show_probability.
function show_probability_Callback(hObject, eventdata, handles)
% hObject    handle to show_probability (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hData
set(handles.axes1);
imagesc(hData.svm.currentPSignal,'Parent',handles.axes1);
colormap('jet')

% --- Executes on button press in show_segmentation.
function show_segmentation_Callback(hObject, eventdata, handles)
% hObject    handle to show_segmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global hData
hData.svm.threshold = str2double(get(handles.seg_threshold,'string'));
hData.svm.currentSegmentation = hData.svm.currentPSignal>hData.svm.threshold;

hData.svm.currentSegmentation = imfill(bwlabel(hData.svm.currentSegmentation));


[CurrentObjNhood,CurrentObjLabels] = bwdist(hData.svm.currentSegmentation );
CurrentObjNhood = double(CurrentObjNhood < 2).*hData.svm.currentSegmentation (CurrentObjLabels);

ed = CurrentObjNhood>0 & hData.svm.currentSegmentation==0; %

R = hData.currentOrigVisualizaion./max(hData.currentOrigVisualizaion(:));
R(edge(hData.currentSegSection))=0;
G=hData.currentOrigVisualizaion./max(hData.currentOrigVisualizaion(:));
G(edge(hData.currentSegSection))=1;
B=hData.currentOrigVisualizaion./max(hData.currentOrigVisualizaion(:));
B(edge(hData.currentSegSection))=1;


%double(hData.svm.currentSegmentation);%
%G = R;%hData.currentOrigVisualizaion./max(hData.currentOrigVisualizaion(:));
%B = R;

R(ed)=1;
G(ed)=0;
B(ed)=0;
RGB = R;
RGB(:,:,2) = G;
RGB(:,:,3) = B;

set(handles.axes1);
imagesc(RGB,'Parent',handles.axes1);
%imagesc(hData.svm.currentSegmentation,'Parent',handles.axes1);
colormap('gray')




% --- Executes on button press in save_clasification.
function save_clasification_Callback(hObject, eventdata, handles)
% hObject    handle to save_clasification (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hData
ProjectName = get(handles.project_name,'String');
hData.ProjectName = ProjectName;
OutputFile = fullfile(hData.PathToSettings,ProjectName);
fprintf('%s: Saving project to\n%s\n',mfilename,OutputFile);
save(OutputFile,'hData','-v7.3');
fprintf('%s: Project saved.\n')

function seg_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to seg_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of seg_threshold as text
%        str2double(get(hObject,'String')) returns contents of seg_threshold as a double


% --- Executes during object creation, after setting all properties.
function seg_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to seg_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function brush_size_Callback(hObject, eventdata, handles)
% hObject    handle to brush_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function brush_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to brush_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function path_settings_Callback(hObject, eventdata, handles)
% hObject    handle to path_settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of path_settings as text
%        str2double(get(hObject,'String')) returns contents of path_settings as a double


% --- Executes during object creation, after setting all properties.
function path_settings_CreateFcn(hObject, eventdata, handles)
% hObject    handle to path_settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in load_settings.
function load_settings_Callback(hObject, eventdata, handles)
% hObject    handle to load_settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%generate a global variable to
global hData

% initialize hData structure
hData = struct();
hData.TraceMouse = false;

% load settings
strPathSettings = get(handles.path_settings,'String');
strProjectName = get(handles.project_name,'String');
structSettings = load(strPathSettings);
[pathstr,name] = fileparts(strPathSettings);

% if hData is found then the setting is an ongoing project
if isfield(structSettings,'hData')
    fprintf('%s: Existing project detected.\n',mfilename)
    hData = structSettings.hData;
    hData.PathToSettings = pathstr;
    hData.SettingsName = name;
    set(handles.project_name,'String',hData.ProjectName);
    set(handles.seg_threshold,'string',num2str(hData.svm.threshold));
    
    
    % update smooth check boxes
    set(handles.SmoothProb,'Value',hData.svm.smooth_p);
    set(handles.SmoothMore,'Value',hData.svm.smooth_pm);
    
    
    set(handles.Section_Number,'String', num2str(hData.currentSectionIndex));
    
    if hData.svm.smooth_p
        set(handles.SmoothMore,'Enable','on');
    else
        set(handles.SmoothMore,'Enable','off');
    end

    
    try
        set(handles.svm_type,'Value',hData.svm.svm_type_val);
    catch
        warning('The SM kernel could not be set. Please set the kernel manually.')
    end
    % if structSettings is found then itis a new project
elseif isfield(structSettings,'structSettings')
    fprintf('%s: New project detected. Loading settings.\n',mfilename)
    % initialize classification variables
    hData.PathToSettings = pathstr;
    hData.SettingsName = name;
    hData.ProjectName = get(handles.project_name,'String');
    
    hData.structSettings = structSettings.structSettings;
    hData.currentImageIndex = 0;
    hData.currentSectionIndex = 1;
    hData.svm.currentClass=false;
    hData.svm.pixeldata=[];
    hData.svm.pixelclass=[];
    hData.svm.currentBackgrLabels = [];
    hData.svm.currentSignalLabels = [];
    hData.svm.currentXTraining=[];
    hData.svm.currentYTraining=[];
    hData.svm.XTraining=[];
    hData.svm.YTraining=[];
else
    warning('Settings file with unknown structure.')
    return
end

hData.svm.currentMethod = 'Eraser';
set(handles.current_color,...
    'BackgroundColor',[0 0 0],...
    'String','Eraser');

fprintf('%s: Loading finished, thanks for waiting.\n',mfilename)

% --- Mouse movement functions.
function moviment_down(hObject, eventdata, handles)
try
    global hData
    hData.TraceMouse=true;
    updateAxes(handles)
catch
end

function movement(hObject, eventdata, handles)
% Unpack gui object
try
    global hData
    if ~hData.TraceMouse
        return
    end
catch
    return
end
updateAxes(handles)

function moviment_up(hObject, eventdata, handles)
%disp('up movement')

try
    global hData
    hData.TraceMouse = false;
    
    
    %xt=round(xt);yt=round(yt);
    indt = hData.svm.currentSignalLabels(:);%sub2ind(size(f1),yt,xt);
    indf = hData.svm.currentBackgrLabels(:);%sub2ind(size(f1),yf,xf);
    
    t=[];f=[];
    for i = 1:length(hData.currentSectionPixelFeatures)
        t = [t hData.currentSectionPixelFeatures{i}(indt)];
        f = [f hData.currentSectionPixelFeatures{i}(indf)];
    end
    
    hData.svm.currentXTraining = [t;f];
    hData.svm.currentYTraining = [ones(size(t,1),1);ones(size(f,1),1)+1];
catch
end

% --- Visualization of current view.
function updateAxes(handles)

global hData

brushSize = round(get(handles.brush_size,'Value')*20);

if brushSize<4
    if mod(brushSize,2) == 0;
        brushSize = brushSize+1;
    end
    SE = strel('square', brushSize);
elseif brushSize<8
    brushSize = floor(brushSize/2);
    SE = strel('diamond', brushSize);
else
    brushSize = floor(brushSize/2);
    SE = strel('disk', brushSize);
end

im = false( size(hData.currentOrigVisualizaion)  );%A(:,:,round(get(a,'Value'))*3+2);

pos = get(handles.axes1,'CurrentPoint');
%flag_btn=get(hObject,'SelectionType');

[m,n]=size(im);
cm=round(pos(1,2));
cn=round(pos(1,1));

im( max(cm,1):min(cm,m) , max(cn,1):min(cn,n) ) = true;
im( cm, cn ) = true;
im = imdilate(im,SE);

switch hData.svm.currentMethod
    case 'Background'
        hData.svm.currentBackgrLabels(im(:)) = true;
        hData.svm.currentSignalLabels(im(:)) = false;
    case 'Signal'
        hData.svm.currentBackgrLabels(im(:)) = false;
        hData.svm.currentSignalLabels(im(:)) = true;
    case 'Eraser'
        hData.svm.currentBackgrLabels(im(:)) = false;
        hData.svm.currentSignalLabels(im(:)) = false;
end

visualizeCurrentImage(handles)

function visualizeCurrentImage(handles)

global hData

ed = edge(hData.currentSegSection);
R = hData.currentOrigVisualizaion./max(hData.currentOrigVisualizaion(:));
R(ed)=0;
R(hData.svm.currentBackgrLabels(:))=max(R(:));
R(hData.svm.currentSignalLabels(:))=0;

G = hData.currentOrigVisualizaion./max(hData.currentOrigVisualizaion(:));
G(ed)=1;
G(hData.svm.currentSignalLabels(:))=max(G(:));
G(hData.svm.currentBackgrLabels(:))=0;

B = hData.currentOrigVisualizaion./max(hData.currentOrigVisualizaion(:));
B(ed)=1;
B(hData.svm.currentBackgrLabels(:))=0;
B(hData.svm.currentSignalLabels(:))=0;

hData.currentRGB(:,:,1) = R;
hData.currentRGB(:,:,2) = G;
hData.currentRGB(:,:,3) = B;
hData.currentRGB=hData.currentRGB./max(hData.currentRGB(:));

set(handles.axes1);
imagesc(hData.currentRGB,'Parent',handles.axes1);

function project_name_Callback(hObject, eventdata, handles)
% hObject    handle to project_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of project_name as text
%        str2double(get(hObject,'String')) returns contents of project_name as a double


% --- Executes during object creation, after setting all properties.
function project_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to project_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in svm_type.
function svm_type_Callback(hObject, eventdata, handles)
% hObject    handle to svm_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns svm_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from svm_type


% --- Executes during object creation, after setting all properties.
function svm_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to svm_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function im = rescale_image(im,lb,ub)

lb = quantile(im(:),lb);
ub = quantile(im(:),ub);

im(im>ub)=ub;
im = im-lb;
im(im<0)=0;


% --------------------------------------------------------------------
function Close_Callback(hObject, eventdata, handles)
% hObject    handle to Close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in show_orig_im_segmentation.
function show_orig_im_segmentation_Callback(hObject, eventdata, handles)
% hObject    handle to show_orig_im_segmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of show_orig_im_segmentation



function im = addCellNucSeg(im)

global hData
[root_path,imName] = fileparts(hData.structSettings.pathToImages{hData.currentImageIndex});
imName = [imName '.png'];
[root_path] = fileparts(root_path);
seg_path = fullfile(root_path,'SEGMENTATION');

%replace
warning('Asuming Action number is 1, and channel is 1 for nuclei and cells segmentation.')
expression = 'A(\d+)Z';
segName = regexprep(imName,expression,'A01Z','preservecase');
expression = 'C(\d+).png';
nucSegName = regexprep(segName,expression,'C01_SegmentedNuclei.png','preservecase');
celSegName = regexprep(segName,expression,'C01_SegmentedCells.png','preservecase');

R=im(:,:,1);
G=im(:,:,2);
B=im(:,:,3);

try
    imSegNuc = edge(imread(fullfile(seg_path,nucSegName)));
    R(imSegNuc>0)=0;
    G(imSegNuc>0)=1;
    B(imSegNuc>0)=0;
catch
    warning('No nuclei segmentation found.')
end

try
    imSegCel = edge(imread(fullfile(seg_path,celSegName)));
    R(imSegCel>0)=1;
    G(imSegCel>0)=0;
    B(imSegCel>0)=1;
catch
    warning('No cell segmentation found.')
end


im(:,:,1)=R;
im(:,:,2)=G;
im(:,:,3)=B;



function LowQuantile_Callback(hObject, eventdata, handles)
% hObject    handle to LowQuantile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LowQuantile as text
%        str2double(get(hObject,'String')) returns contents of LowQuantile as a double


% --- Executes during object creation, after setting all properties.
function LowQuantile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LowQuantile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function HighQuantile_Callback(hObject, eventdata, handles)
% hObject    handle to HighQuantile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of HighQuantile as text
%        str2double(get(hObject,'String')) returns contents of HighQuantile as a double


% --- Executes during object creation, after setting all properties.
function HighQuantile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HighQuantile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SmoothProb.
function SmoothProb_Callback(hObject, eventdata, handles)
% hObject    handle to SmoothProb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hData.svm.smooth_p = get(handles.SmoothProb,'Value');
if hData.svm.smooth_p
    set(handles.SmoothMore,'Enable','on');
else
    set(handles.SmoothMore,'Enable','off');
end
% Hint: get(hObject,'Value') returns toggle state of SmoothProb


% --- Executes on button press in refresh.
function refresh_Callback(hObject, eventdata, handles)
% hObject    handle to refresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hData

hData.lb = str2double( get(handles.LowQuantile,'String') );
hData.ub = str2double( get(handles.HighQuantile,'String') );
lb = quantile(hData.currentOrigImage(:),hData.lb);
ub = quantile(hData.currentOrigImage(:),hData.ub);
%hData.ub = quantile(hData.currentOrigImage(:),0.99995);
hData.currentOrigVisualizaion = hData.currentOrigSection-lb;
hData.currentOrigVisualizaion(hData.currentOrigVisualizaion<0)=0;
hData.currentOrigVisualizaion(hData.currentOrigVisualizaion>ub)=ub;

visualizeCurrentImage(handles)


% --- Executes on button press in SmoothMore.
function SmoothMore_Callback(hObject, eventdata, handles)
% hObject    handle to SmoothMore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SmoothMore



function Section_Number_Callback(hObject, eventdata, handles)
% hObject    handle to Section_Number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Section_Number as text
%        str2double(get(hObject,'String')) returns contents of Section_Number as a double


% --- Executes during object creation, after setting all properties.
function Section_Number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Section_Number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Optimize_Weights.
function Optimize_Weights_Callback(hObject, eventdata, handles)
% hObject    handle to Optimize_Weights (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

warning('Optimization may take several hours. Press "Crtl+c" to cancel')
pause(5)

global hData
XTraining = hData.svm.XTraining;
YTraining = hData.svm.YTraining;
kernel = hData.svm.kernel; 
TestLabel = rand(length(YTraining),1)>0.5;

lb = ones(1,size(XTraining,2)) .* 0.01;
ub = ones(1,size(XTraining,2)) .* 6;

if ~isfield(hData.svm,'best_weights')
    hData.svm.best_weights = ones(1,size(XTraining,2));
end

if size(hData.svm.best_weights,1)>1
    hData.svm.best_weights = ones(1,size(XTraining,2));
end

InitialPopulation = repmat(hData.svm.best_weights,50,1).*(rand(50,length(hData.svm.best_weights))+0.5);
InitialPopulation(1,:)=hData.svm.best_weights;

opts = gaoptimset('PopulationSize', 50, 'Generations', 50,'InitialPopulation',InitialPopulation);

fun = @(x) pcsf.svm_optimize_weights(x,XTraining,YTraining,TestLabel,kernel);

[xbest, fbest, exitflag] = ga(fun, size(XTraining,2), [], [], [], [], ...
    lb, ub, [], [], opts);

hData.svm.best_weights = xbest;


% --- Executes on button press in PlotFeatures.
function PlotFeatures_Callback(hObject, eventdata, handles)
% hObject    handle to PlotFeatures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global hData
for i=1:length(hData.currentSectionPixelFeatures)
    figure;
    imagesc(hData.currentSectionPixelFeatures{i})
    colorbar
    title(sprintf('Feature %d.4',i))
end
