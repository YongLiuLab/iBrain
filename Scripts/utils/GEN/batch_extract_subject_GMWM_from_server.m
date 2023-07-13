clear all;clc;
dataRoot=[filesep,'data',filesep,'yhzhang',filesep,'T1',filesep,'mri',filesep];
BNatlas_struct=load_nii('rBN_Atlas_246_1mm.nii');
BNatlas_data=BNatlas_struct.img;
ROI_num=length(unique(BNatlas_data))-1;
load('ADNI1_orig.mat');
subject_num=length(ADNI1.subname);
all_subject_Atlas_WM=zeros(subject_num,ROI_num);
all_subject_WM=zeros(subject_num,1);
all_subject_GM=zeros(subject_num,1);
for temp_subject=1:subject_num   
    disp(['Current processing subject: ', ADNI1.subname{temp_subject},32, num2str(temp_subject),' out of ', num2str(subject_num)])
    temp_InputStruct=dir([dataRoot,'mwp2*',ADNI1.subname{temp_subject},'*.nii']);
    temp_subject_WM_path=[dataRoot,temp_InputStruct(1).name];
    temp_subject_WM_struct=load_nii(temp_subject_WM_path);
    temp_subject_WM_data=temp_subject_WM_struct.img;
    for temp_ROI=1:ROI_num
        temp_ROI_WMs = temp_subject_WM_data(find(BNatlas_data==temp_ROI));
        all_subject_Atlas_WM(temp_subject,temp_ROI)=sum(temp_ROI_WMs(:));
    end
    all_subject_WM(temp_subject)=sum(sum(sum(temp_subject_WM_data)));

    temp_InputStruct=dir([dataRoot,'mwp1*',ADNI1.subname{temp_subject},'*.nii']);
    temp_subject_GM_path=[dataRoot,temp_InputStruct(1).name];
    temp_subject_GM_struct=load_nii(temp_subject_GM_path);
    temp_subject_GM_data=temp_subject_GM_struct.img;
    all_subject_GM(temp_subject)=sum(sum(sum(temp_subject_GM_data)));
end
ADNI1.Atlas_WM=all_subject_Atlas_WM;
ADNI1.All_WM=all_subject_WM;
ADNI1.All_GM=all_subject_GM;
save('ADNI1_orig.mat','ADNI1','-append');
disp('All subject process finished...')