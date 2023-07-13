function fea=jt_main_feature_test(t2imglist,maskimglist,out_dir)
% pause(7200)
% warning off %#ok<WNOFF>
resdir = strcat(out_dir);
if isempty(dir(resdir))
    mkdir(resdir)
end
% t2imglist = textread(input_txt,'%s');
% maskimglist =  textread(input_mask,'%s');
fea_nm_1st ={'energy';'energy1';'kurtosis';'maximum';'mean';'mad';'median';'minimum';'range';...
    'rms';'skewness';'std';'entropy';'uniformity';'var'};%15 features
fea_nm_2nd ={'Area';'Volume';'Compactness1';'Compactness2';'Max3dDiam';'SphericalDisprop';'Spherity';'Surf2VolRatio';};%8 features

fea_nm_3rd = { 'Autocorrelation';'ClusterProminence';'ClusterShade';'ClusterTendency';...
    'Contrast';'Correlation';'DifferenceEntropy';'Dissimilarity';'Energy';'Entropy';'Homogeneity1';'Homogeneity2';'IMC1';'IMC2';...
    'IDMN';'IDN';'InverseVariance';'MaximumProbability';'SumAverage';'SumEntropy';'SumVariance';'Variance';'ShortRunEmphasis';...
    'LongRunEmphasis';'GrayLevelNonuniformity';'RunLengthNonuniformity';'RunPercentage';'LowGrayLevelRunEmphasis';'HighGrayLevelRunEmphasis';...
    'ShortRunLowGrayLevelEmphasis';'ShortRunHighGrayLevelEmphasis';'LongRunLowGrayLevelEmphasis';'LongRunHighGrayLevelEmphasis';...
    };%33 features



for i = 1:length(t2imglist)
    i
    [pth t2nm]= (fileparts(t2imglist{i}));
    [pth masknm] = (fileparts(maskimglist{i}));
    if ~strcmp(t2nm,masknm)
        fprintf('please double check the subject %d with name %s, for the t2 img and mask img not match\n',i,t2nm);
    end
    t2_nm = t2imglist{i};
    mask_nm = maskimglist{i};
    
    [tu_img t2_img roi_img] = jt_read_mask(mask_nm,t2_nm);
    t2_img = t2_img/max(t2_img(:));
    [x_minax,y_minax,z_minax] = cut_img(tu_img);
    tu_img = tu_img(x_minax(1):x_minax(2),y_minax(1):y_minax(2),z_minax(1):z_minax(2));
    t2_img=t2_img(x_minax(1):x_minax(2),y_minax(1):y_minax(2),z_minax(1):z_minax(2));
    fea.f1(i) = jt_1st_feature(t2_img,tu_img);
    fea.f2(i) =jt_2nd_ShapeSize(t2_img,tu_img);
    fea.f3(i)=jt_3rd_Texture(t2_img,tu_img,32,16);
end

