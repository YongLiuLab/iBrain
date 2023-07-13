function plot_Yeo7_bar(risk_in_Yeo7_system,save_name)
colors = [228 26 28;55 126 184;77 175 74;152 78 163;255 127 0;166 86 40;247 129 191] ./ 255;  

figure;
h = bar(risk_in_Yeo7_system);
h.FaceColor = 'flat';
for i = 1:numel(risk_in_Yeo7_system)
    h.CData(i,:) = colors(i,:);
%     set(h(i),'FaceColor',colors(i,:));
end

Yeo7_names={'Visual','Smatomotor','Dorsal Attention','Ventral Attention','Limbic','Frontal Parietal','Default Mode'};
set(gca, 'XTick', 1:numel(risk_in_Yeo7_system), 'XTickLabels', Yeo7_names);
saveas(gcf,save_name);
close(gcf);
end