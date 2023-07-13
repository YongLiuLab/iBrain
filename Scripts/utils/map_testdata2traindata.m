function norm_data = map_testdata2traindata(test_data,train_map,flag)
switch flag 
    case 1
        for i = 1:length(test_data(:,1,1))
            norm_data(i,:,:) = (squeeze(test_data(i,:,:))-squeeze(train_map.R2SN_m(:,:,2)))./(train_map.R2SN_m(:,:,1)-train_map.R2SN_m(:,:,2));
        end
    case 2
        for i = 1:length(test_data(:,1,1))
            norm_data(i,:,:) = (squeeze(test_data(i,:,:))-squeeze(train_map.RM_m(:,:,2)))./(train_map.RM_m(:,:,1)-train_map.RM_m(:,:,2));
        end
    case 3
        for i = 1:length(test_data(:,1))
            norm_data(i,:) = (test_data(i,:)'-train_map.GM_m(:,2))./(train_map.GM_m(:,1)-train_map.GM_m(:,2));
        end
    case 4
       for i = 1:length(test_data(:,1))
            norm_data(i,:) = (test_data(i,:)'-train_map.MR2SN_m(:,2))./(train_map.MR2SN_m(:,1)-train_map.MR2SN_m(:,2));
       end 
    
end