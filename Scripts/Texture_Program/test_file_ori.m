%
clear;clc
addpath('D:\Dr_work\subtype_validation\Re1203\Code\Texture_Program')
atlas_path  = 'D:\Dr_work\subtype_validation\Re1203\atlas_choose\';
atlas_name = 'rSchaefer2018_100Parcels_7Networks_order_FSLMNI152_1mm.nii';
atlas_img = load_nii([atlas_path,atlas_name]);
max_atlas_num = max(atlas_img.img(:));
img_list = textread('D:\Dr_work\subtype_validation\Re1203\Code\Texture_plus\test_file\test_go\img.txt','%s');
lab_list = [];
for i = 1:length(img_list)
    lab_list = [lab_list;{['D:\Dr_work\subtype_validation\Re1203\Code\temp\temp.nii']}];
end
for i = 1
    z1 = zeros(size(atlas_img.img));
    z1(find(atlas_img.img==i))=1;
    if sum(z1(:))>500
        nn = atlas_img;
        nn.img=z1;
        save_nii(nn,'D:\Dr_work\subtype_validation\Re1203\Code\temp\temp.nii');
        out_dir = ['D:\Dr_work\subtype_validation\Re1203\Code\Texture_plus\test_file\test_go\out_ori\',atlas_name(1:end-4),'\roi_',sprintf('%03d',i),'\'];
        jt_main_feature(img_list,lab_list,out_dir)
    end
end