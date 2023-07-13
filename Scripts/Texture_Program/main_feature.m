clear;clc
out_path = 'E:\M_center_usrBN\Fea_out\MRP_S\';
img_in = 'E:\M_center_usrBN\Fea_out\MRP_S\img\img.txt';
mask_in = 'E:\M_center_usrBN\label\';
img_list = textread('D:\Dr_work\R2SN_pleline\pyradiomics_test\img\img.txt','%s');
for i = 1:93
    lab  = {};
    for j = 1:length(img_list)
        lab = [lab;{['D:\Dr_work\R2SN_pleline\pyradiomics_test\label\roi',sprintf('%03d',i),'.nii']}];
    end
    jt_main_feature(img_list,lab,['D:\Dr_work\R2SN_pleline\pyradiomics_test\out_mat\roi_',sprintf('%03d',i)])
end