clear;clc
addpath(genpath('/data/liugroup/home/kzhao/software/brant-stable'));
out_path = '/data/liugroup/home/kzhao/R2SN_pleline/out/';
img_in = '/data/liugroup/home/kzhao/R2SN_pleline/ATNs_norm/img.txt';
img_list = textread(img_in,'%s');
for i = 50:93
    lab  = {};
    for j = 1:length(img_list)
        lab = [lab;{['/data/liugroup/home/kzhao/R2SN_pleline/label/roi',sprintf('%03d',i),'.nii']}];
    end
    jt_main_feature(img_list,lab,[out_path,'roi_',sprintf('%03d',i)])
end