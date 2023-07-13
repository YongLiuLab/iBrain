function [TP,TN,FN,FP] = conclusion(score,test_label,wi)
pre(find(score>wi))=1;
pre(find(score<=wi))=0;
TP = length(intersect(find(pre==0),find(test_label==0)));
TN = length(intersect(find(pre==1),find(test_label==1)));
FN = length(intersect(find(pre==0),find(test_label==1)));
FP = length(intersect(find(pre==1),find(test_label==0)));