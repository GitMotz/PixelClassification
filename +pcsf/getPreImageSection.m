function imSection = getPreImageSection(im,PadSize,sectionNumber)


tim = padarray(im,[PadSize PadSize],'symmetric','both');

imSize=size(im);
sectionSize = round(imSize./5);
sectionIndex = zeros(5,5);

sectionIndex(:)=1:length(sectionIndex(:));
[i,j]=find(sectionIndex==sectionNumber);

YMIN = sectionSize(1)*(i-1)+1;
XMIN = sectionSize(2)*(j-1)+1;
RECT = [XMIN YMIN sectionSize(2)+PadSize*2 sectionSize(1)+PadSize*2];

imSection = imcrop(tim,RECT);
