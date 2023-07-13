function []= jt_reslice(old_fn,new_fn,voxelsize,method)
% for read img and reslice using reslice_nii
% method (optional)  -	1, 2, or 3
%  			1:  for Trilinear interpolation
%  			2:  for Nearest Neighbor interpolation
%  			3:  for Fischer's Bresenham interpolation
%  			'method' is 1 if it is default or empty.
% for more help read reslice_nii by NIFTI data format can be found on:
% http://nifti.nimh.nih.gov     - Jimmy Shen (jshen@research.baycrest.org)

% Yong Liu, yliu@nlpr.ia.ac.cn; 

if nargin <2
    new_fn = strcat(old_fn(1:end-3),'_res.nii');
end
if nargin <3
   voxelsize = [1 1 1];
end
if nargin <4
    method =1;
end
reslice_nii(old_fn,new_fn,voxelsize,1,[],method);