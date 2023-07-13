function generate_python_T1_report_input(train_data,test_data,GM_img,atlas_img,hippo_reference_GM_path,atlas_reference_GM_path,IsReportToIndSpace,orig_img_path,img_save_name)
ibrainScore = ibrain_computing(train_data,test_data,true);
savePath=fileparts(img_save_name);
save([savePath,filesep,'iBrainScore.mat'],'ibrainScore');

iBrainPath = fileparts(which('iBrain.m'));
load([iBrainPath,filesep,'model_data',filesep,'BN246_Yeo7_map_indexs.mat']);
Yeo_names={'Sub_Cortical','Visual','Smatomotor','Dorsal_Attention','Ventral_Attention','Limbic','Frontal_Parietal','Default_Mode'};

load([iBrainPath,filesep,'model_data',filesep,'ADNI1.mat']);
NC_indexes=find(strcmp(ADNI1.group,'CN'));
NC_Atlas_GM = ADNI1.GM(NC_indexes,:);
mean_NC_Atlas_GM=mean(NC_Atlas_GM);
std_NC_Atlas_GM=std(NC_Atlas_GM);
test_atlas_GM=(test_data.GM-mean_NC_Atlas_GM)./std_NC_Atlas_GM;%convert to mean 0 std1 distribution
Risk_thres=-1.96;%95% confidence interval
risk_test_atlas_feature_indexs = find(test_atlas_GM<Risk_thres);%only focus on significant decreaed volumes,mean-2sigma
%disp(['Subject risk brain region indexs: ',num2str(risk_test_atlas_feature_indexs)])

opts = detectImportOptions(atlas_reference_GM_path);
opts.Encoding = 'GB2312';
Atlas_reference_data = readtable(atlas_reference_GM_path, opts);
all_belongings_names= Atlas_reference_data.feature;
opts = detectImportOptions(hippo_reference_GM_path);
opts.Encoding = 'GB2312';
hippo_reference_data = readtable(hippo_reference_GM_path, opts);
all_hippo_belongings_names =hippo_reference_data.feature;

saveStruct=struct();
saveStruct.ibrainScore=ibrainScore;
%generate hippocampus feature
hippo_ROI_indexs=216:219;
for temp_feature=1:length(hippo_ROI_indexs)
    eval(['saveStruct.',all_hippo_belongings_names{temp_feature},'=test_data.GM(hippo_ROI_indexs(temp_feature));']);
end

%generate risk features by grey matter volume decrease
if ~isempty(risk_test_atlas_feature_indexs)
    all_risk_feature_name=cell(length(risk_test_atlas_feature_indexs),1);
    %generate risk image brain map and top risk image volume features for python report
    output_img=GM_img;
    output_img_data=zeros(size(atlas_img.img));
    for temp_ROI=1:length(risk_test_atlas_feature_indexs)
        output_img_data(find(atlas_img.img==risk_test_atlas_feature_indexs(temp_ROI)))=test_atlas_GM(risk_test_atlas_feature_indexs(temp_ROI));
        all_risk_feature_name{temp_ROI}=all_belongings_names{risk_test_atlas_feature_indexs(temp_ROI)};
    end
    %save risk brain image data
    output_img.img=output_img_data;
    %output_img.fileprefix=img_save_name;
    if exist(img_save_name,'file')
        if isunix
            system(['rm -f ',32, img_save_name]);
        else
            system(['del /f /q ',32,img_save_name]);
        end
    end    
    save_nii(output_img,img_save_name,0);
    %add coregistration to individual space if required
    if IsReportToIndSpace
        [filepath, ~, ~] = fileparts(orig_img_path);
        system(['python ',[iBrainPath,filesep,'Scripts',filesep,'ANTsCoreg.py'],32, '--ref_img=',...
            orig_img_path,32,'--input_img=',...
            [iBrainPath,filesep,'Template',filesep,'cat12_MNI152_T1_1mm.nii'],32,'--input_atlas_mask=',...
            img_save_name,32,'--savePath=',filepath]);          
    end

    for temp_risk_feature_index=1:length(risk_test_atlas_feature_indexs)
        eval(['saveStruct.',all_risk_feature_name{temp_risk_feature_index},'=test_data.GM(risk_test_atlas_feature_indexs(temp_risk_feature_index));']);
    end
    %generate Yeo 7 network mapping alteration statistics, offer optentional instructions for TDCS  
    temp_subject_allcluster_map_to_Yeo_num=zeros(7,1);
    for temp_ROI=1:length(temp_subject_allcluster_map_to_Yeo_num)
        temp_subject_allcluster_map_to_Yeo_num(temp_ROI)=length(find(ROI_Yeo7_belongings(risk_test_atlas_feature_indexs)==temp_ROI));
    end
    savePath=fileparts(img_save_name);
    bar_Save_name=[savePath,filesep,'T1_risk_map_to_Yeo7.jpg'];
    plot_Yeo7_bar(temp_subject_allcluster_map_to_Yeo_num,bar_Save_name);

    %correlate and FDR correct correlation values for risk map and meta
    %maps, generate heatmap based on correlation results.
    heat_save_name=[savePath,filesep,'T1_risk_corr_with_cognition.jpg'];
    risk_map=zeros(length(test_atlas_GM),1);
    risk_map(risk_test_atlas_feature_indexs)=test_data.GM(risk_test_atlas_feature_indexs)';
    generate_heatmap_for_risk_and_neurosynth(risk_map,heat_save_name);   
end

savetable=struct2table(saveStruct);
savePath=fileparts(img_save_name);
writetable(savetable,[savePath,filesep,'features.csv']);

ROI_num=max(unique(atlas_img.img));
AddsaveStruct=struct();
%generate whole brain volume features
for temp_ROI=1:ROI_num
    eval(['AddsaveStruct.',all_belongings_names{temp_ROI},'=test_data.GM(temp_ROI);']);
end
Addsavetable=struct2table(AddsaveStruct);
writetable(Addsavetable,[savePath,filesep,'features_wholebrain.csv']);


