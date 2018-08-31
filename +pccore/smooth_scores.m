function [Score_1,Score_2] = smooth_scores(hData,Score_1,Score_2)

if isfield(hData.svm,'smooth_p_param')
    h = fspecial('gaussian', hData.svm.smooth_p_param(1), hData.svm.smooth_p_param(2)); h = h./sum(h(:));
else
    warning('No parameters for smoothing found. Assuming high level of smoothing.')
    h = fspecial('gaussian', 5, 1.5); h = h./sum(h(:));
end

Score_1 = imfilter(Score_1,h,'symmetric');
Score_2 = imfilter(Score_2,h,'symmetric');

switch hData.svm.kernel
    case 'GentleBoost'
    otherwise
        Score_1=Score_1./(Score_1+Score_2);
        Score_2=Score_2./(Score_1+Score_2);
end
