%
clear;clc
addpath('/data/kzhao/kzhao/Subtype_AS/Re_1203/Code/Texture_Program')
addpath(genpath('/data/kzhao/kzhao/Subtype_AS/Re_1203/Code/Code/brant-master2'))
atlas_path  = '/data/kzhao/kzhao/Subtype_AS/Re_1203/Atlas/';
atlas_name = 'rBN_Atlas_246_1mm.nii';
atlas_img = load_nii([atlas_path,atlas_name]);
max_atlas_num = max(atlas_img.img(:));
img_list = textread('/data/kzhao/Neuroimage_T1/MCAD/Reg/AD_S01.txt','%s');
lab_list = [];
for i = 1:length(img_list)
    lab_list = [lab_list;{['/data/kzhao/kzhao/Subtype_AS/Re_1203/Code/Code/temp/temp.nii']}];
end
for i = 1:246
    z1 = zeros(size(atlas_img.img));
    z1(find(atlas_img.img==i))=1;
    nn = atlas_img;
    nn.img=z1;
    save_nii(nn,'/data/kzhao/kzhao/Subtype_AS/Re_1203/Code/Code/temp/temp.nii');
    out_dir = ['/data/kzhao/NAT_Test/ibrain/MCADI/AD_S01'];
    R_fea(i)=jt_main_feature_test(img_list,lab_list,out_dir);
end
save(out_dir,'R_fea')