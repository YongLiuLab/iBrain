iBrainPath=fileparts(which('iBrain.m'));
Atlas_Path=[iBrainPath,filesep,'Atlas',filesep,'BN_246_Atlas_1mm.nii'];
Atlas_struct=load_nii(Atlas_Path);
Atlas_data=Atlas_struct.img;
%hippo_indexs=[216,217,218,219];%in BN246 atlas
hippo_indexs=477:482;%in Whole512 atlas
hippo_data=zeros(size(Atlas_data));
for temp_ROI_index=1:length(hippo_indexs)
    hippo_data(find(Atlas_data==hippo_indexs(temp_ROI_index)))=1;
end
output_struct=Atlas_struct;
save_name=[iBrainPath,filesep,'Atlas',filesep,'hippo_mask.nii'];
output_struct.img=hippo_data;
save_nii(output_struct,save_name,0);


