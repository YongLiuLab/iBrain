function ibrainscore = ibrain_computing(train_data,test_data,norm_test_flag)
%add tolerance for NaN values
train_data.R2SN(find(isnan(train_data.R2SN)))=0;
train_data.MR2SN(find(isnan(train_data.MR2SN)))=0;
train_data.RF(find(isnan(train_data.RF)))=0;
test_data.R2SN(find(isnan(test_data.R2SN)))=0;
test_data.MR2SN(find(isnan(test_data.MR2SN)))=0;
test_data.RF(find(isnan(test_data.RF)))=0;

iBrainPath = fileparts(which('iBrain.m'));
load([iBrainPath,filesep,'model_data',filesep,'indices.mat']);
load([iBrainPath,filesep,'model_data',filesep,'parm.mat']);
load([iBrainPath,filesep,'model_data',filesep,'tr_score.mat']);
load([iBrainPath,filesep,'model_data',filesep,'acc_inner.mat']);
load([iBrainPath,filesep,'model_data',filesep,'train_map.mat']);

if nargin<3
    norm_test_flag=true;
end

for fea_index = 1:4
    if norm_test_flag
        switch fea_index
            case 1
                train_data1 = train_data.R2SN;train_label = train_data.label;
                test_data1 = map_testdata2traindata(test_data.R2SN,train_map,1);test_label = 1;
            case 2
                train_data1 = train_data.RF;train_label = train_data.label;
                test_data1 = map_testdata2traindata(test_data.RF,train_map,2);test_label = 1;
            case 3
                train_data1 = train_data.GM;train_label = train_data.label;
                test_data1 = map_testdata2traindata(test_data.GM,train_map,3);test_label = 1;
            case 4
                train_data1 = train_data.MR2SN;train_label = train_data.label;
                test_data1 = map_testdata2traindata(test_data.MR2SN,train_map,4);test_label = 1;
        end
    else
        switch fea_index
            case 1
                train_data1 = train_data.R2SN;train_label = train_data.label;
                test_data1 = test_data.R2SN;test_label = 1;
            case 2
                train_data1 = train_data.RF;train_label = train_data.label;
                test_data1 = test_data.RF;test_label = 1;
            case 3
                train_data1 = train_data.GM;train_label = train_data.label;
                test_data1 = test_data.GM;test_label = 1;
            case 4
                train_data1 = train_data.MR2SN;train_label = train_data.label;
                test_data1 = test_data.MR2SN;test_label = 1;
        end
    end

    if fea_index<=2
        for multi_indx = 1:size(train_map.GM_m,1)
            c4 = squeeze(train_data1(:,multi_indx,:));
            c4a = squeeze(test_data1(:,multi_indx,:))';
            c4a(find(isnan(c4a)))=0;
            bestc = parm(fea_index,multi_indx,1);bestg = parm(fea_index,multi_indx,2);
            model = svmtrain(train_label,c4,[' -c  ',num2str(2^bestc),'  -g  ',num2str(2^bestg),' -q  ']);
            [predict_label,accuracy, decision_value] = svmpredict(test_label,c4a,model);
            testscore_inner(:,multi_indx)=decision_value;
        end
        trainscore_inner = squeeze(tr_score(:,fea_index,:));
        acc_pre_1 = squeeze(acc_inner(fea_index,:));
        [bestacc,bestc,bestg,bestscore]=svm_cgchoose_r2sn(trainscore_inner(:,find(acc_pre_1>prctile(acc_pre_1,75))),train_label,indices);
        model = svmtrain(train_label,trainscore_inner(:,find(acc_pre_1>prctile(acc_pre_1,75)) ),[' -c  ',num2str(2^bestc),'  -g  ',num2str(2^bestg),' -q  ']);
        [predict_label,accuracy, decision_value] = svmpredict(test_label,testscore_inner(:,find(acc_pre_1>prctile(acc_pre_1,75))),model);
        trainscore(:,fea_index)=bestscore;testscore(:,fea_index)=decision_value;
    else
        [bestacc,bestc,bestg,bestscore]=svm_cgchoose_r2sn(train_data1,train_label,indices);
        model = svmtrain(train_label,train_data1,[' -c  ',num2str(2^bestc),'  -g  ',num2str(2^bestg),' -q  ']); 
        test_data1(find(isnan(test_data1)))=0;
        [predict_label,accuracy, decision_value] = svmpredict(test_label,double(test_data1),model);
        trainscore(:,fea_index)=bestscore;testscore(:,fea_index)=decision_value;
    end
end
theta = glmfit(trainscore(:,[1:4]), train_label, 'binomial', 'link', 'logit');
% ibrainscore = glmval(theta,testscore(:,[1:4]),'logit');
ibrainscore = theta(1)+testscore(:,1:4)*theta(2:5);