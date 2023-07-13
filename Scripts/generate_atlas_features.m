function RF = generate_atlas_features(inputFilePath,atlasPath)
atlas_Data=load_nii(atlasPath);
ROI_num = max(atlas_Data.img(:));
RF=zeros(ROI_num,47);%47 represents sum of 1st order and 3rd order feature number
for ROI_index=1:ROI_num
    temp_Mask_data = zeros(size(atlas_Data.img));
    temp_Mask_data(find(atlas_Data.img==ROI_index))=1;
    [x_minax,y_minax,z_minax] = cut_img(temp_Mask_data);
    temp_Mask_data = temp_Mask_data(x_minax(1):x_minax(2),y_minax(1):y_minax(2),z_minax(1):z_minax(2));

    temp_img= load_nii(inputFilePath);
    temp_img_data = temp_img.img;
    temp_img_data = temp_img_data/max(temp_img_data(:));
    temp_img_data = temp_img_data(x_minax(1):x_minax(2),y_minax(1):y_minax(2),z_minax(1):z_minax(2));
    features_1st = jt_1st_feature(temp_img_data,temp_Mask_data);
    features_1st = struct2cell(features_1st);
    features_3rd = jt_3rd_Texture(temp_img_data,temp_Mask_data,32,16);
    features_3rd = struct2cell(features_3rd);
    idx = cellfun(@(x) any(isnan(x)) || isempty(x), features_1st);
    features_1st(idx)={single(0)};
    features_1st = cellfun(@(x) single(x), features_1st);
    idx = cellfun(@(x) any(isnan(x)) || isempty(x), features_3rd);
    features_3rd(idx)={single(0)}; 
    features_3rd = cellfun(@(x) single(x), features_3rd);
    RF(ROI_index,:)=[features_1st;features_3rd];
end
end

