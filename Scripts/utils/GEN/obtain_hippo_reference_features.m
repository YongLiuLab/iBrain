iBrainPath = fileparts(which('iBrain.m'));
load([iBrainPath,filesep,'model_data',filesep,'ADNI1.mat']);
ad_subject_indexs=find(strcmp(ADNI1.group,'Dementia')==1);
nc_subject_indexs=find(strcmp(ADNI1.group,'CN')==1);
hippo_Atlas_indexs=216:219;

hippo_ad_volumes=ADNI1.GM(ad_subject_indexs,hippo_Atlas_indexs);
mean_hippo_ad_volumes=squeeze(mean(hippo_ad_volumes));
std_hippo_ad_volumes=squeeze(std(hippo_ad_volumes));
hippo_ad_volumes_range={};
all_hippo_ad_volumes=sum(hippo_ad_volumes,2);
mean_all_hippo_ad_volumes=mean(all_hippo_ad_volumes);
std_all_hippo_ad_volumes=std(all_hippo_ad_volumes);
all_hippo_ad_volumes_range=[num2str(mean_all_hippo_ad_volumes-std_all_hippo_ad_volumes),'-',num2str(mean_all_hippo_ad_volumes+std_all_hippo_ad_volumes)];

hippo_nc_volumes=ADNI1.GM(nc_subject_indexs,hippo_Atlas_indexs);
mean_hippo_nc_volumes=squeeze(mean(hippo_nc_volumes));
std_hippo_nc_volumes=squeeze(std(hippo_nc_volumes));
hippo_nc_volumes_range={};
all_hippo_nc_volumes=sum(hippo_nc_volumes,2);
mean_all_hippo_nc_volumes=mean(all_hippo_nc_volumes);
std_all_hippo_nc_volumes=std(all_hippo_nc_volumes);
all_hippo_nc_volumes_range=[num2str(mean_all_hippo_nc_volumes-std_all_hippo_nc_volumes),'-',num2str(mean_all_hippo_nc_volumes+std_all_hippo_nc_volumes)];

for temp_ROI=1:length(hippo_Atlas_indexs)
    hippo_ad_volumes_range=[hippo_ad_volumes_range;[num2str(mean_hippo_ad_volumes(temp_ROI)-std_hippo_ad_volumes(temp_ROI)),'-',...
        num2str(mean_hippo_ad_volumes(temp_ROI)+std_hippo_ad_volumes(temp_ROI))]];
    hippo_nc_volumes_range=[hippo_nc_volumes_range;[num2str(mean_hippo_nc_volumes(temp_ROI)-std_hippo_nc_volumes(temp_ROI)),'-',...
        num2str(mean_hippo_nc_volumes(temp_ROI)+std_hippo_nc_volumes(temp_ROI))]];
end

mean_hippo_nc_volumes=[mean_hippo_nc_volumes mean_all_hippo_nc_volumes];
std_hippo_nc_volumes=[std_hippo_nc_volumes std_all_hippo_nc_volumes];
hippo_nc_volumes_range=[hippo_nc_volumes_range;all_hippo_nc_volumes_range];
mean_hippo_ad_volumes=[mean_hippo_ad_volumes mean_all_hippo_ad_volumes];
std_hippo_ad_volumes=[std_hippo_ad_volumes std_all_hippo_ad_volumes];
hippo_ad_volumes_range=[hippo_ad_volumes_range;all_hippo_ad_volumes_range];

hippo_volume_feature_names={'Hipp_L_2_1','Hipp_R_2_1','Hipp_L_2_2','Hipp_R_2_2','Hipp_All'};
hippo_volume_feature_chinese_names={'左侧海马1区体积','右侧海马1区体积','左侧海马2区体积','右侧海马2区体积','海马整体体积'};
variable_names={'feature','feature_chinese_name','NC-mean','NC-std','NC-range','AD-mean','AD-std','AD-range'};
hippo_volume_table=table(hippo_volume_feature_names',hippo_volume_feature_chinese_names',mean_hippo_nc_volumes',std_hippo_nc_volumes',hippo_nc_volumes_range,...
    mean_hippo_ad_volumes',std_hippo_ad_volumes',hippo_ad_volumes_range,'VariableNames',variable_names);
writetable(hippo_volume_table,[iBrainPath,filesep,'Report_template',filesep,'hippo_BN246_volume_reference.csv'],'Encoding','GB2312');