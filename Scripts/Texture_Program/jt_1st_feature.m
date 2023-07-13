% compute the 1st order feature 
% refer:  
% Aerts HJ, Velazquez ER, Leijenaar RT, et al. Decoding tumour 
% phenotype by noninvasive imaging using a quantitative radiomics approach. Nature communications 2014;5:4006.
function fea = jt_1st_feature(t2img,maskimg)


if size(t2img) ~= size(maskimg)
    error('please check your mask img and t2 img, for they do not have the same size..\n');
else
   
    roi_img = maskimg.* t2img;
end

%%%% check the mask size also be feature volume size
temp = sum(maskimg(:));
if temp<=50
    fprintf('mask size is too small\n');
end

%%  %% 1st order feature base on 
% fea.vol = temp;
roi_vec = roi_img(find(roi_img));
% fea.size = length(roi_vec);
fea.energy = sum(roi_vec.^2)/length(roi_vec);
fea.kurtosis = kurtosis(roi_vec);
fea.maximum = max(roi_vec);
fea.mean = mean(roi_vec);
fea.mad = mad(roi_vec);
fea.median = median(roi_vec);
fea.minimum = min(roi_vec);
fea.range = range(roi_vec);
fea.rms = rms(roi_vec);
fea.skewness = skewness(roi_vec);
fea.std = std(roi_vec);
fea.var = var(roi_vec);
[fea.entropy,fea.uniformity] = jt_entropy(roi_vec,100);
fea.entropy=single(fea.entropy);
fea.uniformity=single(fea.uniformity);