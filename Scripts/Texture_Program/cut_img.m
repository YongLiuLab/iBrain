%
function [x_minax,y_minax,z_minax] = cut_img(img)
% input_img = load_nii('D:\Dr_work\subtype_validation\Re1203\Code\Texture_plus\test_file\roi_001.nii');
% img = input_img.img;
wi_z = sum(sum(img));
indx = find(wi_z);
z_minax=[indx(1)-1,indx(end)+1];

wi_y = sum(sum(img,3));
indx = find(wi_y);
y_minax=[indx(1)-1,indx(end)+1];

wi_x = sum(sum(img,3),2);
indx = find(wi_x);
x_minax=[indx(1)-1,indx(end)+1];



% nii = make_nii(img_plus);
% save_nii(nii,'test.nii')


