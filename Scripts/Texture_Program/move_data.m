re_path = 'F:\BN\MPR_pre\';

for i = 1:41
    X=sprintf('%03d',i);
    dir_na = strcat(re_path,'Vs',X,'\','b_s',X,'\c1reg_b_s',X,'.nii');
    out = strcat('G:\MPR_S\img\s_',X,'.nii');
    copyfile(dir_na,out)
  


end