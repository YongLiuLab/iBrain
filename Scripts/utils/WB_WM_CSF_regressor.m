function timecourse = WB_WM_CSF_regressor(data)
iBrainPath=fileparts(which('iBrain.m'));
L=size(data,4);
a=spm_vol(strcat(iBrainPath,filesep,'Template',filesep,'mask',filesep,'BrainMask_05_61x73x61.img'));
b=spm_read_vols(a);
wb= b(:)>0;
a=spm_vol(strcat(iBrainPath,filesep,'Template',filesep,'mask',filesep,'WhiteMask_09_61x73x61.img'));
b=spm_read_vols(a);
wm= b(:)>0;
a=spm_vol(strcat(iBrainPath,filesep,'Template',filesep,'mask',filesep,'CsfMask_07_61x73x61.img'));
b=spm_read_vols(a);
csf= b(:)>0;
for f=1:L
    temp=data(:,:,:,f);
    timecourse(f,1)=mean(temp(wb));
    timecourse(f,2)=mean(temp(wm));
    timecourse(f,3)=mean(temp(csf));
end
