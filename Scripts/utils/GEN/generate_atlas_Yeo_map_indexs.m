iBrainPath=fileparts(which('iBrain.m'));
%load Atlas
Atlas_struct=load_nii([iBrainPath,filesep,'Atlas',filesep,'BN_Atlas_246_3mm.nii']);
ROI_num=length(unique(Atlas_struct.img))-1;%remove 0
%load Yeo 7 Atlas
Yeo7_struct=load_nii([iBrainPath,filesep,'Atlas',filesep,'Yeo2011_7Networks_MNI152_FreeSurferConformed3mm.nii']);
%For each ROI in Atlas, define its Yeo 7 belonging by standard of maximum
%overlap with Yeo 7 System
ROI_Yeo7_belongings=zeros(ROI_num,1);
for temp_ROI=1:ROI_num
    Yeo7_max_overlap=0;
    Yeo7_max_overlap_index=0;
    temp_ROI_indexs=find(Atlas_struct.img==temp_ROI);
    for temp_map_ROI=1:7
       temp_Yeo7_indexs=find(Yeo7_struct.img==temp_map_ROI);
       temp_max_overlap=length(intersect(temp_ROI_indexs,temp_Yeo7_indexs));
       if temp_max_overlap>Yeo7_max_overlap
           Yeo7_max_overlap=temp_max_overlap;
           Yeo7_max_overlap_index=temp_map_ROI;
       end    
    end
    ROI_Yeo7_belongings(temp_ROI)=Yeo7_max_overlap_index;%0 represents subcortical
end
%Save result indexs to model_data folder
save([iBrainPath,filesep,'model_data',filesep,'BN246_Yeo7_map_indexs.mat']);