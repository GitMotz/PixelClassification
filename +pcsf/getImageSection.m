function imSection = getImageSection(im,sectionNumber)

imSize=size(im);
sectionSize = round(imSize./5);
sectionIndex = zeros(5,5);
sectionIndex(:)=1:length(sectionIndex(:));
[i,j]=find(sectionIndex==sectionNumber);

YMIN = sectionSize(1)*(i-1)+1;
XMIN = sectionSize(2)*(j-1)+1;
RECT = [XMIN YMIN sectionSize(2) sectionSize(1)];

imSection = imcrop(im,RECT);
