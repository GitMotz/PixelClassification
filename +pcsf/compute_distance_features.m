function [d1,d2,d3,norm_int] = compute_distance_features(orig_image,seg_image)

orig_seg_im = seg_image;

S = sort(orig_image(:),'descend');
orig_image(orig_image>S(11))=S(11);


se = strel('disk',6);
seg_image = imerode(seg_image,se);
ids = unique(seg_image(:));ids(ids==0)=[];
d2 = zeros(size(orig_image));
d3 = zeros(size(orig_image));
norm_int = zeros(size(orig_image));


h = fspecial('gaussian',5,1);
d2_orig =  imfilter(orig_image,h);


h = fspecial('gaussian',10, 3); h = h./sum(h(:));
d3_orig = deconvlucy(d2_orig, h);


d1 = bwdist(seg_image)./20;






for i = 1:length(ids)
    f = seg_image==ids(i);
    im = d2_orig(f(:));
    imt = d3_orig(f(:));
    
    
    im = im-min(im);
    immax = max(im);
    im=im./immax;
    
    im = im .* 256;
    im2 = zeros(size(orig_image));
    im2(f(:))=im;
    level = graythresh(round(im)).*256 .* 1.2; %quantile(im,0.90);
    
    d2(im2>level & f) = 1;
    
    
    imt = imt-min(imt);
    immaxt = max(imt);
    imt=imt./immaxt;
    
    imt = imt .* 256;
    imt2 = zeros(size(orig_image));
    imt2(f(:))=imt;
    levelt = quantile(imt,0.93);%graythresh(round(imt)).*256;
    
    d3(imt2>levelt  & f) = 1;
    
    
    f = orig_seg_im == ids(i);
    
    pix = d2_orig(f(:));
    q1 = quantile(pix,0.005);
    pix = pix-q1;
    q2 = quantile(pix,0.995);
    pix = pix./q2;
    pix(pix>1)=1;
    pix(pix<0)=0;
    
    norm_int(f(:))=pix;
    
end

d2 = bwlabel(d2);


areas = regionprops(d2,'Area');
areas =cat(1,areas(:).Area);
f = find(areas<12);
f = ismember(d2(:),f);
d2(f)=0;
%figure;imagesc(d2)
d2 = bwdist(d2)./20;

d3 = bwlabel(d3);
areas = regionprops(d3,'Area');
areas =cat(1,areas(:).Area);
f = find(areas<16);
f = ismember(d3(:),f);
d3(f)=0;
%figure;imagesc(d2)
d3 = bwdist(d3)./20;



