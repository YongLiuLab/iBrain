clear all;clc;
%% ADNI
dataRoot=[filesep,'data',filesep,'yhzhang',filesep,'T1',filesep,'ADNI_T1',filesep];
Atlas_struct=load_nii('rWhole512_1mm.nii');
Atlas_data=Atlas_struct.img;
ROI_num=length(unique(Atlas_data))-1;
load('ADNI1_orig.mat');
all_subject_indexs=find(strcmp(ADNI1.group,'CN') | strcmp(ADNI1.group,'Dementia'));
subject_num=length(all_subject_indexs);
all_subject_Atlas_WM=zeros(subject_num,ROI_num);
all_subject_Atlas_GM=zeros(subject_num,ROI_num);
all_subject_WM=zeros(subject_num,1);
all_subject_GM=zeros(subject_num,1);

load('wi_90.mat');
template_Data='MNI152_T1_1mm.nii';
R2SN=zeros(subject_num,ROI_num,ROI_num);
MR2SN=zeros(subject_num,ROI_num);
GM=zeros(subject_num,ROI_num);
RF=zeros(subject_num,ROI_num,47);
label=zeros(subject_num,1);
for temp_subject=1:subject_num   
    disp(['Current processing subject: ', ADNI1.subname{all_subject_indexs(temp_subject)},32, num2str(temp_subject),' out of ', num2str(subject_num)])
    %% construct ADNI1_orig saving data
    temp_InputStruct=dir([dataRoot,'mri',filesep,'mwp2*',ADNI1.subname{all_subject_indexs(temp_subject)},'*.nii']);
    temp_subject_WM_path=[dataRoot,'mri',filesep,temp_InputStruct(1).name];
    temp_subject_WM_struct=load_nii(temp_subject_WM_path);
    temp_subject_WM_data=temp_subject_WM_struct.img;
    for temp_ROI=1:ROI_num
        temp_ROI_WMs = temp_subject_WM_data(find(Atlas_data==temp_ROI));
        all_subject_Atlas_WM(temp_subject,temp_ROI)=sum(temp_ROI_WMs(:));
    end
    all_subject_WM(temp_subject)=sum(sum(sum(temp_subject_WM_data)));

    temp_InputStruct=dir([dataRoot,'mri',filesep,'mwp1*',ADNI1.subname{all_subject_indexs(temp_subject)},'*.nii']);
    temp_subject_GM_path=[dataRoot,'mri',filesep,temp_InputStruct(1).name];
    temp_subject_GM_struct=load_nii(temp_subject_GM_path);
    temp_subject_GM_data=temp_subject_GM_struct.img;
    for temp_ROI=1:ROI_num
        temp_ROI_GMs = temp_subject_GM_data(find(Atlas_data==temp_ROI));
        all_subject_Atlas_GM(temp_subject,temp_ROI)=sum(temp_ROI_GMs(:));
    end
    all_subject_GM(temp_subject)=sum(sum(sum(temp_subject_GM_data)));
    %% construct train_data saving data
    temp_InputStruct=dir([dataRoot,'DNR',filesep,'*',ADNI1.subname{all_subject_indexs(temp_subject)},'*.nii']);
    temp_InputFile=[dataRoot,'DNR',filesep,temp_InputStruct(1).name];
    RF(temp_subject,:,:) = generate_atlas_features(temp_InputFile,'Whole512_1mm.nii'); 
    GM_volume = get_GM(temp_subject_GM_struct,Atlas_struct);
    test_data = constructR2SN(squeeze(RF(temp_subject,:,:)),GM_volume,wi); 
    R2SN(temp_subject,:,:)=test_data.R2SN;
    MR2SN(temp_subject,:)=test_data.MR2SN;
    GM(temp_subject,:)=test_data.GM;   
    if strcmp(ADNI1.group(all_subject_indexs(temp_subject)),'CN')
        label(temp_subject)=1;
    end
end

R2SN(find(isnan(R2SN)))=0;
MR2SN(find(isnan(MR2SN)))=0;
RF(find(isnan(RF)))=0;

ADNI1.WM=all_subject_Atlas_WM;
ADNI1.GM=all_subject_Atlas_GM;
ADNI1.All_WM=all_subject_WM;
ADNI1.All_GM=all_subject_GM;
ADNI1.net_plus=R2SN;
ADNI1.RM = RF;
save('ADNI1_Whole512.mat','ADNI1');

train_data.RF=data_norm(RF);
train_data.GM=mapminmax(GM',0,1)';
train_data.R2SN=data_norm(R2SN);
train_data.MR2SN=mapminmax(MR2SN',0,1)';
save('train_data_Whole512.mat','train_data');

%% MCAD
% total_center=8;
% dataRoot=[filesep,'data',filesep,'yhzhang',filesep,'T1',filesep,'MCAD_T1',filesep];
% Atlas_struct=load_nii('rWhole512_1mm.nii');
% Atlas_data=Atlas_struct.img;
% ROI_num=length(unique(Atlas_data))-1;
% load('MCAD_orig.mat');
% subject_num=length(AIBL.group);
% all_subject_Atlas_WM=zeros(subject_num,ROI_num);
% all_subject_Atlas_GM=zeros(subject_num,ROI_num);
% all_subject_WM=zeros(subject_num,1);
% all_subject_GM=zeros(subject_num,1);
% 
% load('wi_90.mat');
% template_Data='MNI152_T1_1mm.nii';
% R2SN=zeros(subject_num,ROI_num,ROI_num);
% MR2SN=zeros(subject_num,ROI_num);
% GM=zeros(subject_num,ROI_num);
% RF=zeros(subject_num,ROI_num,47);
% 
% temp_subject_in_all_center_index=0;
% for temp_center=1:total_center
%     temp_center_subject_indexs=find(AIBL.center==temp_center);
%     for temp_subject=1:length(temp_center_subject_indexs)
%         temp_subject_in_all_center_index=temp_subject_in_all_center_index+1;
%         disp(['Current processing subject: ', AIBL.subname{temp_center_subject_indexs(temp_subject)},32, num2str(temp_subject_in_all_center_index),' out of ', num2str(subject_num)])
%         %% construct ADNI1_orig saving data
%         temp_InputStruct=dir([dataRoot,'AD_S0',num2str(temp_center),filesep,'mri',filesep,'mwp2*',AIBL.subname{temp_center_subject_indexs(temp_subject)},'*.nii']);
%         temp_subject_WM_path=[dataRoot,'AD_S0',num2str(temp_center),filesep,'mri',filesep,temp_InputStruct(1).name];
%         temp_subject_WM_struct=load_nii(temp_subject_WM_path);
%         temp_subject_WM_data=temp_subject_WM_struct.img;
%         for temp_ROI=1:ROI_num
%             temp_ROI_WMs = temp_subject_WM_data(find(Atlas_data==temp_ROI));
%             all_subject_Atlas_WM(temp_subject_in_all_center_index,temp_ROI)=sum(temp_ROI_WMs(:));
%         end
%         all_subject_WM(temp_subject_in_all_center_index)=sum(sum(sum(temp_subject_WM_data)));
% 
%         temp_InputStruct=dir([dataRoot,'AD_S0',num2str(temp_center),filesep,'mri',filesep,'mwp1*',AIBL.subname{temp_center_subject_indexs(temp_subject)},'*.nii']);
%         temp_subject_GM_path=[dataRoot,'AD_S0',num2str(temp_center),filesep,'mri',filesep,temp_InputStruct(1).name];
%         temp_subject_GM_struct=load_nii(temp_subject_GM_path);
%         temp_subject_GM_data=temp_subject_GM_struct.img;
%         for temp_ROI=1:ROI_num
%             temp_ROI_GMs = temp_subject_GM_data(find(Atlas_data==temp_ROI));
%             all_subject_Atlas_GM(temp_subject_in_all_center_index,temp_ROI)=sum(temp_ROI_GMs(:));
%         end
%         all_subject_GM(temp_subject_in_all_center_index)=sum(sum(sum(temp_subject_GM_data)));
%         %% construct train_data saving data
%         temp_InputStruct=dir([dataRoot,'DNR',filesep,'AD_S0',num2str(temp_center),filesep,'*',AIBL.subname{temp_center_subject_indexs(temp_subject)},'*.nii']);
%         temp_InputFile=[dataRoot,'DNR',filesep,'AD_S0',num2str(temp_center),filesep,temp_InputStruct(1).name];
%         RF(temp_subject_in_all_center_index,:,:) = generate_atlas_features(temp_InputFile,'Whole512_1mm.nii');
%         GM_volume = get_GM(temp_subject_GM_struct,Atlas_struct);
%         test_data = constructR2SN(squeeze(RF(temp_subject_in_all_center_index,:,:)),GM_volume,wi);
%         R2SN(temp_subject_in_all_center_index,:,:)=test_data.R2SN;
%         MR2SN(temp_subject_in_all_center_index,:)=test_data.MR2SN;
%         GM(temp_subject_in_all_center_index,:)=test_data.GM;
%     end
% end
% 
% 
% AIBL.WM=all_subject_Atlas_WM;
% AIBL.GM=all_subject_Atlas_GM;
% AIBL.All_WM=all_subject_WM;
% AIBL.All_GM=all_subject_GM;
% AIBL.RM=RF;
% AIBL.net_plus=R2SN;
% AIBL.MR2SN=MR2SN;
% save('AIBL_Whole512.mat','AIBL');


disp('All subject process finished...')