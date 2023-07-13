function gm = get_GM(data,Atlas)
ROI_num=max(unique(Atlas.img));
gm=zeros(1,ROI_num);
data=data.img;
for temp_atlas_region = 1:ROI_num
    gm_k = data(find(Atlas.img==temp_atlas_region));
    gm(temp_atlas_region)=sum(gm_k(:));
end
