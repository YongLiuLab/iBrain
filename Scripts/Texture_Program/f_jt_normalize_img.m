% function f_jt_normalize_img()

warning off %#ok<WNOFF>
t2imglist = textread('I:\MRI_Glioma_Tiantan\Tumor_mask_201601\t2list.txt','%s');
maskimglist =  textread('I:\MRI_Glioma_Tiantan\Tumor_mask_201601\masklist.txt','%s');

mask_all = 0;
for i = 1:length(t2imglist)
    i
    [pth t2nm]= (fileparts(fileparts(t2imglist{i})));
    [pth masknm] = (fileparts(fileparts(maskimglist{i})));
    if ~strcmp(t2nm,masknm)
        fprintf('please double check the subject %d with name %s, for the t2 img and mask img not match\n',i,t2nm);
    end
    t2_nm = t2imglist{i};
    mask_nm = maskimglist{i};
    
    [mask_img t2_img roi_img] = jt_read_mask(mask_nm,t2_nm);
    
%     mask_all = mask_all+mask_img;

    Temp = t2_img - roi_img;
    
    I = find(Temp>=80);
    Temp = Temp(I);
    
   [H Bin] = hist(Temp,500,'nodisplay');
   
   [I J] = sort(H);
   
   Peak_img(i) = mean(Bin(J(end-2:end)));
   
  t2_nm = t2imglist{i};

   [nii]  = load_untouch_nii(t2_nm);
   [pth nm] = fileparts(t2_nm);
   
   nii.fileprefix  = strcat(pth,'\nor_',nm,'.nii');
   nii.img = nii.img.*(1000/Peak_img(i));
   
    save_untouch_nii(nii,nii.fileprefix);
   
 
    
   
end
% mask_nm = maskimglist{1};

[nii]  = load_untouch_nii(mask_nm);
% nii.fileprefix = 'mask_all.nii';
% nii.img = mask_all;
% save_nii(nii,nii.fileprefix,0);
% 
% 
% options = statset('Display','final');
% obj = fitgmdist(H,2);
