function test_data = constructR2SN(temp_load,GM_volume,model_weight)
flag_sub = mapminmax(temp_load.var')';
R2SN(1,:,:)=corr(flag_sub(:, find(model_weight))');
R2SN(find(isnan(R2SN)))=0;%correct NaN value for small ROIs
test_data.GM = GM_volume;
test_data.R2SN = R2SN;
test_data.RF(1,:,:) = temp_load.var;
test_data.MR2SN=squeeze(mean(R2SN))';
end

