function [data_no,t_data_m]=data_norm(tr_data)
ROI_num=size(tr_data,2);
for i = 1:ROI_num
    t_data = squeeze(tr_data(:,i,:));
    t_data_m(i,:,1) = max(t_data);
    t_data_m(i,:,2) = min(t_data);
    data_no(:,i,:)=mapminmax(t_data',0,1)';
end