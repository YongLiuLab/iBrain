function [bestacc,bestc,bestg,bestscore]=svm_cgchoose_r2sn(data,label,sites_tr)
% grid search 
bestacc = 0;
bestc = 0;
bestg = 0;
bestscore=zeros(length(label),1);
indices = unique(sites_tr);
for c = -3:5
    for g = -5:-1
        for isite = 1:length(indices)
            train_data = data(find(sites_tr~=indices(isite)),:);
            test_data = data(find(sites_tr==indices(isite)),:);
            train_label = label(find(sites_tr~=indices(isite)));
            test_label = label(find(sites_tr==indices(isite)));
            model = svmtrain(train_label,train_data,[' -c  ',num2str(2^c),'  -g  ',num2str(2^g),' -q  ']);
            [predict_label,accuracy, decision_value] = svmpredict (test_label,test_data,model);
            acc_cg(isite) = accuracy(1);
            score(find(sites_tr==indices(isite)))=decision_value;
        end
        acc_cg = acc_cg(find(acc_cg));
        if mean(acc_cg)>bestacc
            bestscore = score;
            bestc = c;
            bestg = g;
            bestacc = mean(acc_cg);
        end
    end
end

