function [mask_img t2_img roi_img] = jt_read_mask(maskname,t2name)
% for read mask img and get the t2 img

% maskname = 'M:\works\work_JT\testData\ltu.nii';
% t2name = 'M:\works\work_JT\testData\t2.nii';
[maskimg]  = load_nii(maskname);
[t2img]  = load_nii(t2name);

if maskimg.hdr.dime.dim ~= t2img.hdr.dime.dim
    error('please check your mask img and t2 img, for they do not have the same size..\n');
else
    mask_img = double(maskimg.img>=0.8);
    t2_img = double(t2img.img);
    roi_img = mask_img.* t2_img;
end

% xyz= [maskimg.hdr.hist.srow_x;maskimg.hdr.hist.srow_y;maskimg.hdr.hist.srow_y];



