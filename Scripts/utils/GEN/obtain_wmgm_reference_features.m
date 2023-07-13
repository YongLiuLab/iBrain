iBrainPath = fileparts(which('iBrain.m'));
load([iBrainPath,filesep,'model_data',filesep,'ADNI1.mat']);

ad_subject_indexs=find(strcmp(ADNI1.group,'Dementia')==1);
nc_subject_indexs=find(strcmp(ADNI1.group,'CN')==1);


GM_ad_volumes=ADNI1.All_GM(ad_subject_indexs,:);
mean_GM_ad_volumes=squeeze(mean(GM_ad_volumes));
std_GM_ad_volumes=squeeze(std(GM_ad_volumes));
GM_ad_volumes_range=[num2str(mean_GM_ad_volumes-std_GM_ad_volumes),'-',num2str(mean_GM_ad_volumes+std_GM_ad_volumes)];

GM_nc_volumes=ADNI1.All_GM(nc_subject_indexs,:);
mean_GM_nc_volumes=squeeze(mean(GM_nc_volumes));
std_GM_nc_volumes=squeeze(std(GM_nc_volumes));
GM_nc_volumes_range=[num2str(mean_GM_nc_volumes-std_GM_nc_volumes),'-',num2str(mean_GM_nc_volumes+std_GM_nc_volumes)];

WM_ad_volumes=ADNI1.All_WM(ad_subject_indexs,:);
mean_WM_ad_volumes=squeeze(mean(WM_ad_volumes));
std_WM_ad_volumes=squeeze(std(WM_ad_volumes));
WM_ad_volumes_range=[num2str(mean_WM_ad_volumes-std_WM_ad_volumes),'-',num2str(mean_WM_ad_volumes+std_WM_ad_volumes)];

WM_nc_volumes=ADNI1.All_WM(nc_subject_indexs,:);
mean_WM_nc_volumes=squeeze(mean(WM_nc_volumes));
std_WM_nc_volumes=squeeze(std(WM_nc_volumes));
WM_nc_volumes_range=[num2str(mean_WM_nc_volumes-std_WM_nc_volumes),'-',num2str(mean_WM_nc_volumes+std_WM_nc_volumes)];

GMper_ad_volumes=ADNI1.All_GM(ad_subject_indexs,:)./(ADNI1.All_GM(ad_subject_indexs,:)+ADNI1.All_WM(ad_subject_indexs,:));
mean_GMper_ad_volumes=squeeze(mean(GMper_ad_volumes));
std_GMper_ad_volumes=squeeze(std(GMper_ad_volumes));
GMper_ad_volumes_range=[num2str(mean_GMper_ad_volumes-std_GMper_ad_volumes),'-',num2str(mean_GMper_ad_volumes+std_GMper_ad_volumes)];

GMper_nc_volumes=ADNI1.All_GM(nc_subject_indexs,:)./(ADNI1.All_GM(nc_subject_indexs,:)+ADNI1.All_WM(nc_subject_indexs,:));
mean_GMper_nc_volumes=squeeze(mean(GMper_nc_volumes));
std_GMper_nc_volumes=squeeze(std(GMper_nc_volumes));
GMper_nc_volumes_range=[num2str(mean_GMper_nc_volumes-std_GMper_nc_volumes),'-',num2str(mean_GMper_nc_volumes+std_GMper_nc_volumes)];

WMper_ad_volumes=ADNI1.All_WM(ad_subject_indexs,:)./(ADNI1.All_GM(ad_subject_indexs,:)+ADNI1.All_WM(ad_subject_indexs,:));
mean_WMper_ad_volumes=squeeze(mean(WMper_ad_volumes));
std_WMper_ad_volumes=squeeze(std(WMper_ad_volumes));
WMper_ad_volumes_range=[num2str(mean_WMper_ad_volumes-std_WMper_ad_volumes),'-',num2str(mean_WMper_ad_volumes+std_WMper_ad_volumes)];

WMper_nc_volumes=ADNI1.All_WM(nc_subject_indexs,:)./(ADNI1.All_GM(nc_subject_indexs,:)+ADNI1.All_WM(nc_subject_indexs,:));
mean_WMper_nc_volumes=squeeze(mean(WMper_nc_volumes));
std_WMper_nc_volumes=squeeze(std(WMper_nc_volumes));
WMper_nc_volumes_range=[num2str(mean_WMper_nc_volumes-std_WMper_nc_volumes),'-',num2str(mean_WMper_nc_volumes+std_WMper_nc_volumes)];

Total_ad_volumes=ADNI1.All_GM(ad_subject_indexs,:)+ADNI1.All_WM(ad_subject_indexs,:);
mean_Total_ad_volumes=squeeze(mean(Total_ad_volumes));
std_Total_ad_volumes=squeeze(std(Total_ad_volumes));
Total_ad_volumes_range=[num2str(mean_Total_ad_volumes-std_Total_ad_volumes),'-',num2str(mean_Total_ad_volumes+std_Total_ad_volumes)];

Total_nc_volumes=ADNI1.All_GM(nc_subject_indexs,:)+ADNI1.All_WM(nc_subject_indexs,:);
mean_Total_nc_volumes=squeeze(mean(Total_nc_volumes));
std_Total_nc_volumes=squeeze(std(Total_nc_volumes));
Total_nc_volumes_range=[num2str(mean_Total_nc_volumes-std_Total_nc_volumes),'-',num2str(mean_Total_nc_volumes+std_Total_nc_volumes)];

%concatenate all features
all_NC_mean=[mean_WM_nc_volumes;mean_GM_nc_volumes;mean_WMper_nc_volumes;mean_GMper_nc_volumes;mean_Total_nc_volumes];
all_NC_std=[std_WM_nc_volumes;std_GM_nc_volumes;std_WMper_nc_volumes;std_GMper_nc_volumes;std_Total_nc_volumes];
all_NC_range={WM_nc_volumes_range;GM_nc_volumes_range;WMper_nc_volumes_range;GMper_nc_volumes_range;Total_nc_volumes_range};
all_AD_mean=[mean_WM_ad_volumes;mean_GM_ad_volumes;mean_WMper_ad_volumes;mean_GMper_ad_volumes;mean_Total_ad_volumes];
all_AD_std=[std_WM_ad_volumes;std_GM_ad_volumes;std_WMper_ad_volumes;std_GMper_ad_volumes;std_Total_ad_volumes];
all_AD_range={WM_ad_volumes_range;GM_ad_volumes_range;WMper_ad_volumes_range;GMper_ad_volumes_range;Total_ad_volumes_range};
WMGM_volume_feature_names={'白质体积','灰质体积','白质百分比','灰质百分比','全脑体积'};
variable_names={'feature','NC-mean','NC-std','NC-range','AD-mean','AD-std','AD-range'};

WMGM_volume_table=table(WMGM_volume_feature_names',all_NC_mean,all_NC_std,all_NC_range,all_AD_mean,...
    all_AD_std,all_AD_range,'VariableNames',variable_names);
writetable(WMGM_volume_table,[iBrainPath,filesep,'Report_template',filesep,'WMGM_volume_reference.csv'],'Encoding','GB2312');