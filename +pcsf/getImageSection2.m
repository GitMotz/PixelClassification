function imSection = getImageSection2(im,PadSize)
imSection = im(1+PadSize:end-PadSize,1+PadSize:end-PadSize);
