function cellFilteredImages = GeneratePixelFeatures(OrigImage,FeatureList)

fprintf('%s: Calculating features. This will take a while.\n',mfilename)

if nargin<2
    fprintf('Calculating all features.\n')
    FeatureList = ...
        {'log_5_2';'log_5_3';'log_6_3';'log_8_4';'log_10_5';'log_20_10';...
        'gau_20_10';'gau_9_4';'gau_5_2';...
        'diffgau_1';'diffgau_2';...
        'average_3';'average_5';'average_10';'average_20';...
        'edges';'simple_edges';'edge_log_8_4';'weird_log_8_4';...
        'entropy_log_8_4';'entropy_log_5_2';'loc_back_small';'loc_back_large';...
        'morph_opening';'morph_close';'morph_close_small';'deconv_10_7';'deconv_10_3';'deconv_5_3';...
        'std_dif_1';'zscore_1';'watershed_lines';'watershed_lines_smooth';...
        'watershed_lines_smooth_low';'watershed_lines_smooth_very_low';'tophat_12';'orig_tophat_30'};
end


cellFilteredImages={};

S = sort(OrigImage(:),'descend');
OrigImage(OrigImage>S(11))=S(11);

% get the image at the log scale
cellFilteredImages{1} = log10(OrigImage);
fprintf('.')


for iFeat = 1:length(FeatureList)
    FeatureId = FeatureList{iFeat};
    
    
    switch FeatureId
        
        case 'log_5_2'
            % get laplacian of gaussian filters
            h = fspecial('log', 5, 2);
            cellFilteredImages{end+1} = imfilter(cellFilteredImages{1},h,'symmetric');
            fprintf('.')
            
        case 'log_5_3'
            % get laplacian of gaussian filters
            h = fspecial('log', 5, 3);
            cellFilteredImages{end+1} = imfilter(cellFilteredImages{1},h,'symmetric');
            fprintf('.')
            
        case 'log_6_3'
            % get laplacian of gaussian filters
            h = fspecial('log', 5, 3);
            cellFilteredImages{end+1} = imfilter(cellFilteredImages{1},h,'symmetric');
            fprintf('.')
            
            
        case 'log_8_4'
            h = fspecial('log', 8, 4);
            cellFilteredImages{end+1} = imfilter(cellFilteredImages{1},h,'symmetric');
            fprintf('.')
            
        case 'log_10_5'
            h = fspecial('log', 20, 10);
            cellFilteredImages{end+1} = imfilter(cellFilteredImages{1},h,'symmetric');
            fprintf('.')
            
            
        case 'log_20_10'
            h = fspecial('log', 20, 10);
            cellFilteredImages{end+1} = imfilter(cellFilteredImages{1},h,'symmetric');
            fprintf('.')
            
            
            % get gaussian filters
        case 'gau_20_10'
            h = fspecial('gaussian', 20, 10);
            cellFilteredImages{end+1} = imfilter(cellFilteredImages{1},h,'symmetric');
            fprintf('.')
            
        case 'gau_9_4'
            h = fspecial('gaussian', 9, 4);
            cellFilteredImages{end+1} = imfilter(cellFilteredImages{1},h,'symmetric');
            fprintf('.')
            
        case 'gau_5_2'
            h = fspecial('gaussian', 5, 2);
            cellFilteredImages{end+1} = imfilter(cellFilteredImages{1},h,'symmetric');
            fprintf('.')
            
            % get difference of gaussian filters
        case 'diffgau_1'
            h1 = fspecial('gaussian', 15, 8); h1 = h1./max(h1(:));
            h2 = fspecial('gaussian', 15, 3);  h2 = h2./max(h2(:));
            h = (h1-h2); n = h<0; m = min(h(n(:))); if ~isempty(m); h = h+m; end; h = h./sum(h(:));
            cellFilteredImages{end+1} = imfilter(cellFilteredImages{1},h,'symmetric');
            
        case 'diffgau_2'
            h1 = fspecial('gaussian', 15, 9); h1 = h1./max(h1(:));
            h2 = fspecial('gaussian', 15, 5);  h2 = h2./max(h2(:));
            h = (h1-h2); n = h<0; m = min(h(n(:))); if ~isempty(m); h = h+m; end; h = h./sum(h(:));
            cellFilteredImages{end+1} = imfilter(cellFilteredImages{1},h,'symmetric');
            fprintf('.')
            % h1 = fspecial('gaussian', 40, 25); h1 = h1./max(h1(:));
            % h2 = fspecial('gaussian', 40, 10);  h2 = h2./max(h2(:));
            % h = (h1-h2); n = h<0; m = min(h(n(:))); if ~isempty(m); h = h+m; end; h = h./sum(h(:));
            % cellFilteredImages{end+1} = imfilter(cellFilteredImages{1},h,'symmetric');
            
            
            % get average filters
        case 'average_3'
            h = fspecial('average', 3);
            cellFilteredImages{end+1} = imfilter(cellFilteredImages{1},h,'symmetric');
            fprintf('.')
            
        case 'average_5'
            h = fspecial('average', 5);
            cellFilteredImages{end+1} = imfilter(cellFilteredImages{1},h,'symmetric');
            fprintf('.')
            
        case 'average_10'
            h = fspecial('average', 10);
            cellFilteredImages{end+1} = imfilter(cellFilteredImages{1},h,'symmetric');
            fprintf('.')
            
        case 'average_20'
            h = fspecial('average', 20);
            cellFilteredImages{end+1} = imfilter(cellFilteredImages{1},h,'symmetric');
            fprintf('.')
            
            % get edges
        case 'edges'
            h = fspecial('sobel');
            fedge{1} = imfilter(cellFilteredImages{1},h,'symmetric');%13
            fprintf('.')
            fedge{2} = imfilter(cellFilteredImages{1},h','symmetric');%14
            
            h = fspecial('prewitt');
            fedge{3} = imfilter(cellFilteredImages{1},h,'symmetric');%15
            fprintf('.')
            fedge{4} = imfilter(cellFilteredImages{1},h','symmetric');%16
            
            cellFilteredImages{end+1} = fedge{1};
            cellFilteredImages{end+1} = fedge{2};
            cellFilteredImages{end+1} = fedge{3};
            cellFilteredImages{end+1} = fedge{4};
            
            
            
            % get laplacian of gaussian of edges        
            h = fspecial('laplacian', 0.1);
            cellFilteredImages{end+1} = imfilter(fedge{3},h,'symmetric');
            fprintf('.')
            
            cellFilteredImages{end+1} = imfilter(fedge{4},h,'symmetric');
            fprintf('.')
            
            h = fspecial('log', 30);
            cellFilteredImages{end+1} = imfilter(fedge{3},h,'symmetric');
            fprintf('.')
            cellFilteredImages{end+1} = imfilter(fedge{4},h,'symmetric');
            fprintf('.')
            
            
            % get edges of laplacian of gaussian
        case 'simple_edges'
            
            h = fspecial('sobel');
            fedge{1} = imfilter(cellFilteredImages{1},h,'symmetric');%13
            fprintf('.')
            fedge{2} = imfilter(cellFilteredImages{1},h','symmetric');%14
            
            h = fspecial('prewitt');
            fedge{3} = imfilter(cellFilteredImages{1},h,'symmetric');%15
            fprintf('.')
            fedge{4} = imfilter(cellFilteredImages{1},h','symmetric');%16
            
            cellFilteredImages{end+1} = fedge{1};
            cellFilteredImages{end+1} = fedge{2};
            cellFilteredImages{end+1} = fedge{3};
            cellFilteredImages{end+1} = fedge{4};
            
            
        case 'edge_log_8_4'
            h = fspecial('log', 8, 4);
            tim = imfilter(cellFilteredImages{1},h,'symmetric');
            
            h = fspecial('prewitt');
            cellFilteredImages{end+1} = abs(imfilter(abs(tim),h,'symmetric'));
            fprintf('.')
            cellFilteredImages{end+1} = abs(imfilter(abs(tim),h','symmetric'));
            fprintf('.')
            
            % get extra filters
        case 'weird_log_8_4'
            h = fspecial('log', 8, 4);
            tim = imfilter(cellFilteredImages{1},h,'symmetric');
            h = [1 1 1 1 1;1 0 0 0 1;1 0 0 0 1;1 0 0 0 1;1 1 1 1 1];h = h./sum(h(:));
            cellFilteredImages{end+1} = imfilter(tim,h','symmetric');
            fprintf('.')
            
            
            % get entropies
        case 'entropy_log_8_4'
            h = fspecial('log', 8, 4);
            tim = imfilter(cellFilteredImages{1},h,'symmetric');
            se = strel('disk',5);
            cellFilteredImages{end+1} = entropyfilt(abs(tim),se.getnhood);
            fprintf('.')
            
        case 'entropy_log_5_2'
            h = fspecial('log', 5, 2);
            tim = imfilter(cellFilteredImages{1},h,'symmetric');
            se = strel('disk',5);
            cellFilteredImages{end+1} = entropyfilt(abs(tim),se.getnhood);
            fprintf('.')
            
            % get morphological openings
            %se = strel('disk',3);
            
        case 'morph_opening'
            se = strel('disk',5);
            cellFilteredImages{end+1} = imopen(cellFilteredImages{1},se);
            fprintf('.')
            
            %se = strel('disk',5);
            se = strel('disk',8);
            cellFilteredImages{end+1} = imopen(cellFilteredImages{1},se);
            fprintf('.')
            
            se = strel('disk',10);
            cellFilteredImages{end+1} = imopen(cellFilteredImages{1},se);
            fprintf('.')
            
            
            % get local background
        case 'loc_back_small'
            % a = cellFilteredImages{1};
            % a = colfilt(padarray(a,[10 10],nan,'both'),[10 10],'sliding',@(a) quantile(a,0.1));
            % a = a(11:end-10,11:end-10);
            a = imresize(cellFilteredImages{1},0.5);
            a = colfilt(padarray(a,[5 5],nan,'both'),[5 5],'sliding',@(a) quantile(a,0.1));
            a = imresize(a(6:end-5,6:end-5),2);
            cellFilteredImages{end+1} = a;
            
            
        case 'loc_back_large'
            a = imresize(cellFilteredImages{1},0.25);
            a = colfilt(padarray(a,[5 5],nan,'both'),[10 10],'sliding',@(a) quantile(a,0.05));
            a = imresize(a(6:end-5,6:end-5),4);
            cellFilteredImages{end+1} = a;
            fprintf('.')
            
            
        case 'morph_close'
            se = strel('disk',3);
            cellFilteredImages{end+1} = imclose((cellFilteredImages{1}),se);fprintf('.')
            se = strel('disk',5);
            cellFilteredImages{end+1} = imclose((cellFilteredImages{1}),se);fprintf('.')
            se = strel('disk',7);
            cellFilteredImages{end+1} = imclose((cellFilteredImages{1}),se);fprintf('.')
            se = strel('disk',10);
            cellFilteredImages{end+1} = imclose((cellFilteredImages{1}),se);fprintf('.')
            se = strel('disk',15);
            cellFilteredImages{end+1} = imclose((cellFilteredImages{1}),se);fprintf('.')
            se = strel('disk',20);
            cellFilteredImages{end+1} = imclose((cellFilteredImages{1}),se);fprintf('.')
            
        case 'morph_close_small'
            se = strel('disk',3);
            cellFilteredImages{end+1} = imclose((cellFilteredImages{1}),se);fprintf('.')
            se = strel('disk',5);
            cellFilteredImages{end+1} = imclose((cellFilteredImages{1}),se);fprintf('.')
            
            
        case 'deconv_10_7'
            h = fspecial('gaussian', 10, 7);
            cellFilteredImages{end+1} = deconvlucy(cellFilteredImages{1}, h);
            
        case 'deconv_10_3'
            h = fspecial('gaussian', 10, 3);
            cellFilteredImages{end+1} = deconvlucy(cellFilteredImages{1}, h);
            
        case 'deconv_7_2'
            h = fspecial('gaussian', 7, 2);
            cellFilteredImages{end+1} = deconvlucy(cellFilteredImages{1}, h);
            
        case 'deconv_5_3'
            h = fspecial('gaussian', 5, 3); h = h./sum(h(:));
            cellFilteredImages{end+1} = deconvlucy(cellFilteredImages{1}, h);
            
            
        case 'std_dif_1'
            h = fspecial('gaussian', 7, 2);
            tim = deconvlucy(cellFilteredImages{1}, h);
            
            a = (stdfilt(tim, ones(3)));
            b = (stdfilt(tim, ones(9)));
            
            h = fspecial('gaussian',5,2);
            tim = (imfilter(log(abs(b.^2-a.^2)),h,'symmetric'));
            
            se = strel('disk',5);
            cellFilteredImages{end+1} = ( imclose( tim , se ) );
            
            
        case 'zscore_1'
            tim = padarray(cellFilteredImages{1},[20 20],'symmetric','both');
            h1 = fspecial('gaussian',5,2);
            a = stdfilt(tim, ones(3));
            h = fspecial('average',3);
            b = imfilter(tim, h);
            c = (tim-b)./a;
            c(c>3)=3;
            c(c<-3)=-3;
            c = imfilter(c,h1);
            cellFilteredImages{end+1} = c(21:end-20,21:end-20);
            
            
            a = stdfilt(c, ones(5));
            h = fspecial('average',5);
            b = imfilter(c, h);
            cellFilteredImages{end+1} = a(21:end-20,21:end-20);
            cellFilteredImages{end+1} = b(21:end-20,21:end-20);
            
            
        case 'watershed_lines'
            tim = cellFilteredImages{1};
            imres = watershed(tim);
            cellFilteredImages{end+1} = double(imres==0);fprintf('.');
            for i = 1:5
                
                h = fspecial('gaussian',i*3,i);
                h = h./sum(h(:));
                timt = imfilter(tim,h);
                imres(:,:,i+1) = watershed(timt);
            end
            cellFilteredImages{end+1} = nanmean(double(imres==0),3);fprintf('.');
            
            h = fspecial('gaussian',i*3,i);
            h = h./sum(h(:));
            cellFilteredImages{end+1} = imfilter(nanmean(double(imres==0),3),h);fprintf('.');
            
        case 'watershed_lines_smooth'
            
            tim = cellFilteredImages{1};
            imres = watershed(tim);
            for i = 1:5
                
                h = fspecial('gaussian',i*3,i);
                h = h./sum(h(:));
                timt = imfilter(tim,h);
                imres(:,:,i+1) = watershed(timt);
            end
            h = fspecial('gaussian',i*3,i);
            h = h./sum(h(:));
            cellFilteredImages{end+1} = imfilter(nanmean(double(imres==0),3),h);fprintf('.');
            
            
            %             figure;imagesc(a(21:end-20,21:end-20))
            %             figure;imagesc(c(21:end-20,21:end-20))
            %
            %             se = strel('disk',3);
            %
            %             h = fspecial('gaussian', 15, 5);h=abs(h-max(h(:)));h = h./sum(h(:));
            %             figure;imagesc(imfilter(a(21:end-20,21:end-20),h));
            %
            
            
        case 'watershed_lines_smooth_low'
            
            tim = cellFilteredImages{1};
            imres = watershed(tim);
            for i = 1:5
                
                h = fspecial('gaussian',i*3,i);
                h = h./sum(h(:));
                timt = imfilter(tim,h);
                imres(:,:,i+1) = watershed(timt);
            end
            
            h = fspecial('gaussian',5,2);
            h = h./sum(h(:));
            cellFilteredImages{end+1} = imfilter(nanmean(double(imres==0),3),h);fprintf('.');
            
            
        case 'watershed_lines_smooth_very_low'
            
            tim = cellFilteredImages{1};
            imres = watershed(tim);
            for i = 1:5
                
                h = fspecial('gaussian',i*3,i);
                h = h./sum(h(:));
                timt = imfilter(tim,h);
                imres(:,:,i+1) = watershed(timt);
            end
            
            h = fspecial('gaussian',5,1.2);
            h = h./sum(h(:));
            cellFilteredImages{end+1} = imfilter(nanmean(double(imres==0),3),h);fprintf('.');
            
        case 'tophat_12'
            se = strel('disk',12);
            h = fspecial('gaussian',5,1.5);
            tim_orig =  imfilter(cellFilteredImages{1},h);
            cellFilteredImages{end+1} = imtophat(tim_orig,se);
            
        case 'orig_tophat_30'
            se = strel('disk',30);
            h = fspecial('gaussian',5,0.75);
            tim_orig =  imfilter(OrigImage,h);
            cellFilteredImages{end+1} = imtophat(tim_orig,se)./100;
            
            
            
            
        otherwise
            warning('Feature case unknown.')
    end
end


% for i = 1:length(cellFilteredImages)
%    figure(i);imagesc(cellFilteredImages{i});
% end



fprintf('\n%s: End of calculation. Thanks for waiting.\n',mfilename)
