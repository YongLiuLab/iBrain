clear all;clc
ROI_num=512;
%% ADNI
train_map=struct();
iBrainPath=fileparts(which('iBrain.m'));
load([iBrainPath,filesep,'model_data',filesep,'ADNI1.mat']);
load([iBrainPath,filesep,'model_data',filesep,'indices.mat']);
label_adni(find(strcmp(ADNI1.group,'CN')),1)=1;
label_adni(find(strcmp(ADNI1.group,'Dementia')),1)=3;
R2SN_adni = ADNI1.net_plus;
GM_adni = ADNI1.GM;
RM_adni = ADNI1.RM;
clear ADNI1
MR2SN_adni = squeeze(mean(R2SN_adni,3));
MR2SN_adni(find(isnan(MR2SN_adni)))=0;
indx = find((~isnan(label_adni))&(label_adni~=2)&(label_adni)&(label_adni<4));
label_adni = label_adni(indx);
R2SN_adni = R2SN_adni(indx,:,:);
R2SN_adni(find(isnan(R2SN_adni)))=0;
[R2SN_adni,train_map.R2SN_m] = data_norm(R2SN_adni);
R2SN_adni(find(isnan(R2SN_adni)))=0;
GM_adni = double(GM_adni(indx,:,:));
GM_adni(find(isnan(GM_adni)))=0;
train_map.GM_m(:,1) =max(GM_adni);train_map.GM_m(:,2)  =min(GM_adni);
GM_adni = mapminmax(GM_adni',0,1)';
GM_adni(find(isnan(GM_adni)))=0;
MR2SN_adni = MR2SN_adni(indx,:,:);
train_map.MR2SN_m(:,1)  = max(MR2SN_adni);train_map.MR2SN_m(:,2) = min(MR2SN_adni);
MR2SN_adni = mapminmax(MR2SN_adni',0,1)';
MR2SN_adni(find(isnan(MR2SN_adni)))=0;
RM_adni = RM_adni(indx,:,:);
RM_adni(find(isnan(RM_adni)))=0;
[RM_adni,train_map.RM_m] = data_norm(RM_adni);
RM_adni(find(isnan(RM_adni)))=0;
label_adni(find(label_adni==3))=0;
save('train_map_Whole512.mat','train_map');
%% MCADI
load([iBrainPath,filesep,'model_data',filesep,'AIBL.mat']);
label_aibl = AIBL.label;
% AIBL.GM(687,1)=0;
indx = find((~isnan(label_aibl))&(label_aibl~=2)&(AIBL.GM(:,1))&(label_aibl<4)&(label_aibl>0));
label_aibl = label_aibl(indx);
label_aibl(find(label_aibl==3))=0;
R2SN_aibl = AIBL.net_plus(indx,:,:);
R2SN_aibl(find(isnan(R2SN_aibl)))=0;
R2SN_aibl = map_testdata2traindata(R2SN_aibl,train_map,1);
R2SN_aibl(find(isnan(R2SN_aibl)))=0;
MR2SN_aibl = squeeze(mean(R2SN_aibl,3));
MR2SN_aibl(find(isnan(MR2SN_aibl)))=0;
MR2SN_aibl = map_testdata2traindata(MR2SN_aibl,train_map,4);
MR2SN_aibl(find(isnan(MR2SN_aibl)))=0;
GM_aibl = AIBL.GM(indx,:);
GM_aibl(find(isnan(GM_aibl)))=0;
GM_aibl = map_testdata2traindata(GM_aibl,train_map,3);
GM_aibl(find(isnan(GM_aibl)))=0;
RM_aibl = AIBL.RM(indx,:,:);
RM_aibl(find(isnan(RM_aibl)))=0;
RM_aibl = map_testdata2traindata(RM_aibl,train_map,2);
RM_aibl(find(isnan(RM_aibl)))=0;
clear AIBL 
label_aibl(find(label_aibl==3))=0;
%%
site = [zeros(length(label_aibl),1)+1;zeros(length(label_adni),1)+2];
R2SN = [R2SN_aibl;R2SN_adni];
RM = [RM_aibl;RM_adni];
GM = [GM_aibl;GM_adni];
MR2SN = [MR2SN_aibl;MR2SN_adni];
label = [label_aibl;label_adni];
clear R2SN_aibl R2SN_adni  %reduce memory usage
clear GM_aibl GM_adni 
clear RM_aibl RM_adni
%%
train_indx = find(site==2);
test_indx = find(site==1);

for fea_index = 1:4
     switch fea_index 
            case 1 
                train_data = R2SN(train_indx,:,:);train_label = label(train_indx);
                test_data = R2SN(test_indx,:,:);test_label = label(test_indx); 
            case 2
                train_data = RM(train_indx,:,:);train_label = label(train_indx);
                test_data = RM(test_indx,:,:);test_label = label(test_indx); 
            case 3 
                train_data = double(GM(train_indx,:,:));train_label = label(train_indx);
                test_data = double(GM(test_indx,:,:));test_label = label(test_indx); 
            case 4
                train_data = MR2SN(train_indx,:,:);train_label = label(train_indx);
                test_data = MR2SN(test_indx,:,:);test_label = label(test_indx); 
     end
     if fea_index<=2   
          for multi_indx = 1:ROI_num
               c4 = squeeze(train_data(:,multi_indx,:));
               c4a = squeeze(test_data(:,multi_indx,:));
               [bestacc,bestc,bestg,bestscore]=svm_cgchoose_r2sn(c4,train_label,indices);
               parm(fea_index,multi_indx,1)=bestc;
               parm(fea_index,multi_indx,2)=bestg;
               model = svmtrain(train_label,c4,[' -c  ',num2str(2^bestc),'  -g  ',num2str(2^bestg),' -q  ']); 
               [predict_label,accuracy, decision_value] = svmpredict(test_label,c4a,model); 
               trainscore_inner(:,multi_indx)=bestscore;testscore_inner(:,multi_indx)=decision_value;
               acc_pre_1(multi_indx)=bestacc;
               acc_inner(fea_index,multi_indx)=bestacc;
               tr_score(:,fea_index,multi_indx)=bestscore;
          end
         [bestacc,bestc,bestg,bestscore]=svm_cgchoose_r2sn(trainscore_inner(:,find(acc_pre_1>prctile(acc_pre_1,75))),train_label,indices);
          model = svmtrain(train_label,trainscore_inner(:,find(acc_pre_1>prctile(acc_pre_1,75)) ),[' -c  ',num2str(2^bestc),'  -g  ',num2str(2^bestg),' -q  ']); 
         [predict_label,accuracy, decision_value] = svmpredict(test_label,testscore_inner(:,find(acc_pre_1>prctile(acc_pre_1,75))),model); 
         trainscore(:,fea_index)=bestscore;testscore(:,fea_index)=decision_value;
     else 
           [bestacc,bestc,bestg,bestscore]=svm_cgchoose_r2sn(train_data,train_label,indices);
           model = svmtrain(train_label,train_data,[' -c  ',num2str(2^bestc),'  -g  ',num2str(2^bestg),' -q  ']); 
           [predict_label,accuracy, decision_value] = svmpredict(test_label,test_data,model); 
           trainscore(:,fea_index)=bestscore;testscore(:,fea_index)=decision_value;
         
     end
end
indx = intersect((find(test_label<4)),(find(test_label>=0)));
theta = glmfit(trainscore(:,[1:4]), train_label, 'binomial', 'link', 'logit');
score = glmval(theta,testscore(:,[1:4]),'logit');
score_ori = theta(1)+testscore(:,1:4)*theta(2:5);
[TP,TN,FN,FP] = conclusion(score(indx),test_label(indx),0.5);
SEN1 = TP / (TP + FN);
SPE1 = TN / (TN + FP);
ACC1 = (TP + TN) / (TP + FP + FN + TN);
[auc, curve] = rocplot(score(indx),test_label(indx), 1, 0);

save('ROC_data.mat','trainscore','train_label','theta','score','test_label','indx');
save('parm_Whole512.mat','parm');
save('tr_score_Whole512.mat','tr_score');
save('acc_inner_Whole512.mat','acc_inner');



