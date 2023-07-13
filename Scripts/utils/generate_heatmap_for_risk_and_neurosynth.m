function generate_heatmap_for_risk_and_neurosynth(risk_feature_data,savename)
iBrainPath = fileparts(which('iBrain.m'));
load([iBrainPath,filesep,'model_data',filesep,'neurosynth_comptopics.mat']);
corr_data=zeros(32,1);
for temp_comp=1:length(corr_data)
    corr_data(temp_comp)=corr(comp_region(temp_comp,:)',risk_feature_data);
end
corr_data=reshape(corr_data,8,4);
complabels=reshape(complabels,8,4);

figure;
imagesc(corr_data);
colormap('cool'); 

matrixSize = size(corr_data);

for i = 1:matrixSize(1)
    for j = 1:matrixSize(2)
        value = corr_data(i, j);  
        text(j, i, strcat(num2str(value, '%.2f'), 32, complabels{i,j}), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end

colorbar;
axis off;
saveas(gcf,savename);
close(gcf);

end