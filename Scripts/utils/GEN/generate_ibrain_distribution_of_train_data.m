clear all;clc;
delete(gcp('nocreate'));
parpool('local',8);
if exist('gcp.m','file')
    try
        gcp;
    end
elseif parpool('size')==0
    try
        parpool;
    end
end
iBrainPath = fileparts(which('iBrain.m'));
load([iBrainPath,filesep,'model_data',filesep,'train_data.mat']);
train_subject_num=length(train_data.label);
train_subject_ibrain_score=zeros(train_subject_num,1);
norm_test_flag=false;
for temp_train_subject=1:train_subject_num
    disp(['Current processing ',num2str(temp_train_subject),' out of ',num2str(train_subject_num), ' subjects'])
    temp_test_data=struct();
    temp_test_data.R2SN=train_data.R2SN(temp_train_subject,:,:);
    temp_test_data.MR2SN=mean(temp_test_data.R2SN);
    temp_test_data.GM=train_data.GM(temp_train_subject,:);
    temp_test_data.RF=train_data.RF(temp_train_subject,:,:);
    temp_test_data.label=train_data.label(temp_train_subject);
    train_subject_ibrain_score(temp_train_subject) = 1-ibrain_computing(train_data,temp_test_data,norm_test_flag); 
end
%obtain mean and std for NC and AD group
nc_train_subject_ibrain_score=train_subject_ibrain_score(find(train_data.label(1:train_subject_num)==1));
nc_mean_train_subject_ibrain_score=mean(nc_train_subject_ibrain_score);
nc_std_train_subject_ibrain_score=std(nc_train_subject_ibrain_score);
ad_train_subject_ibrain_score=train_subject_ibrain_score(find(train_data.label(1:train_subject_num)==0));
ad_mean_train_subject_ibrain_score=mean(ad_train_subject_ibrain_score);
ad_std_train_subject_ibrain_score=std(ad_train_subject_ibrain_score);  
save([iBrainPath,filesep,'model_data',filesep,'train_data_ibrain.mat'],'nc_mean_train_subject_ibrain_score',...
    'nc_std_train_subject_ibrain_score','ad_mean_train_subject_ibrain_score','ad_std_train_subject_ibrain_score',...
    'nc_train_subject_ibrain_score','ad_train_subject_ibrain_score');