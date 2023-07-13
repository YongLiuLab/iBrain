clear;clc
addpath('/home/kzhao/Subtype_AS/Re_1203/Code/Code/Texture_Program')
addpath(genpath('/home/kzhao/Subtype_AS/Re_1203/Code/Code/brant-master2'))
atlas_path  = '/home/kzhao/Subtype_AS/Re_1203/Atlas/';
atlas_name = 'rSchaefer2018_800Parcels_7Networks_order_FSLMNI152_1mm.nii';
atlas_img = load_nii([atlas_path,atlas_name]);
max_atlas_num = max(atlas_img.img(:));
img_list = textread('/home/kzhao/Subtype_AS/Re_1203/DATA/img.txt','%s');
lab_list = [];
for i = 1:length(img_list)
    lab_list = [lab_list;{['/home/kzhao/Subtype_AS/Re_1203/Code/Code/temp/temp83.nii']}];
end
for i =[403:409,429:600]
    z1 = zeros(size(atlas_img.img));
    z1(find(atlas_img.img==i))=1;
    nn = atlas_img;
    nn.img=z1;
    save_nii(nn,'/home/kzhao/Subtype_AS/Re_1203/Code/Code/temp/temp83.nii');
    out_dir = ['/home/kzhao/Subtype_AS/Re_1203/Feature_out/',atlas_name(1:end-4),'/roi_',sprintf('%03d',i),'/'];
    jt_main_feature(img_list,lab_list,out_dir)
end