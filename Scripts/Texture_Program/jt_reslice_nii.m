% % function jt_reslice_nii
clear all
clc
% % list = textread('I:\MRI_Glioma_Tiantan\Tumor_mask_201601\mask2list.txt','%s');
% % voxelsize = [1 1 1];
% % method = 2;
% % for i = 1:length(list)
% %     i
% %     
% %     old_fn = list{i};
% %     [pth nm] = fileparts(list{i});    
% %     new_fn = strcat(pth,'\res_',nm,'.nii');    
% %     jt_reslice(old_fn,new_fn,voxelsize,method)
% % end
%%
list = textread('I:\MRI_Glioma_Tiantan\Tumor_mask_201601\list.txt','%s');
voxelsize = [1 1 1];
method = 1;
for i = 1:length(list)
    i
    
    old_fn = list{i};
    [pth nm] = fileparts(list{i});    
    new_fn = strcat(pth,'\res_',nm,'.nii');    
    jt_reslice(old_fn,new_fn,voxelsize,method)
end