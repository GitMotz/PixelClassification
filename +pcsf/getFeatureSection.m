function cellImSections = getFeatureSection(cellIms,sectionNumber)

imSize=size(cellIms{1});
sectionSize = round(imSize./5);
sectionIndex = zeros(5,5);
sectionIndex(:)=1:length(sectionIndex(:));
[i,j]=find(sectionIndex==sectionNumber);

YMIN = sectionSize(1)*(i-1)+1;
XMIN = sectionSize(2)*(j-1)+1;
RECT = [XMIN YMIN sectionSize(2) sectionSize(1)];

cellImSections = cell(size(cellIms));
for i = 1:length(cellIms)
cellImSections{i} = imcrop(cellIms{i},RECT);
end