% nii_name = dir(fullfile('G:\ADNI_N3\','*nii'));
nii_name = textread('ADNI_MRI.txt','%s');
for i = 1:length(nii_name)
    A=cell2mat(nii_name(i));
    [x,y] = fileparts(A);
    out = strcat('I:\ADNI_MRI\',y,'.nii');
    copyfile(A,out)
end