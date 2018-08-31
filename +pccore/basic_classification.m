function [Score_1,Score_2] = basic_classification(hData,CurrentData,CellSeg) 

fprintf('%s: Classifying current section. This will take a moment.\n',mfilename);
tic
segmentation=zeros(size(CurrentData{1}));
Score_1=segmentation;
Score_2=segmentation;
%segmentation(:)=1:length(segmentation(:));


CellSegOn = nargin == 3;
if CellSegOn && any(CellSeg(:)>0)
   [CurrentObjNhood,CurrentObjLabels] = bwdist(CellSeg );
   CellSeg = uint16(CurrentObjNhood < 2).*CellSeg(CurrentObjLabels);
   cc = regionprops(CellSeg,'PixelIdxList');
   PixelList = cat(1,cc(:).PixelIdxList);
   ix = randi(1000,length(PixelList),1);
   [~,~,ix]=unique(ix);
   PixelList = arrayfun(@(i) PixelList(ix==i),unique(ix),'uniformoutput',false);
   clear cc
else
   PixelList = segmentation;
   PixelList(:) = 1:length(segmentation(:));   
   PixelList = arrayfun(@(i) PixelList(:,i),1:size(PixelList,2),'uniformoutput',false);
end


% build the data for classification of the curret segment
allData = nan(length(CurrentData{1}(:)),length(CurrentData));
for i = 1:length(CurrentData)
    allData(:,i) = CurrentData{1}(:);
    CurrentData(1)=[];
end

if isfield(hData.svm,'best_weights');
    for i = 1:size(allData,2)
        allData(:,i)=allData(:,i).*hData.svm.best_weights(i);
    end
end

clear CurrentData

allData=single(allData);


n = 1;
switch hData.svm.kernel
    case 'GentleBoost'
        [a, b]  = predict(hData.svm.CompactSVMModel, allData);
        Score_1(:)  = b(:,1);
        Score_2(:)  = b(:,2);
        
    case 'SVMEnsemble'
        
        for j = 1:length(PixelList)
            i = PixelList{j};
            %[~, p]  = predict(hData.svm.CompactSVMModel1, allData(i,:));
            [a, b]  = pccore.svm_ensemble_predict(hData.svm.CompactSVMModel, [allData(i,:)]);
            segmentation(i) = a(:,:,1);
            b = nanmean(b,3);
            Score_1(i)  = b(:,1);
            Score_2(i)  = b(:,2);
            n = 1+n
        end
        
        
    otherwise
        for j = 1:length(PixelList)
            
            i = PixelList{j};
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

if CellSegOn
    Score_1(CellSeg==0)=min(Score_1(:));
    Score_2(CellSeg==0)=max(Score_2(:));
end

toc