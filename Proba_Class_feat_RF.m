function [ new_x ] = Proba_Class_feat_RF( RFmodel,x )
Nview = length(RFmodel);
class_feature1 = []; proba_feature1 = [];
class_feature2 = []; proba_feature2 = [];
for i=1:Nview
    [predict_label,proba] = predict(RFmodel{i}, x{i}');
    predict_label = cell2mat(predict_label);
    predict_label = str2num(predict_label);
    class_feature1 = [class_feature1;predict_label'];
%     proba = proba(:,1);
    proba_feature1 = [proba_feature1;proba'];
end
% new_x = x;
new_x{1} = class_feature1';
% new_x{end} = NormalizeFea(new_x{end},0);
% new_x{end} = new_x{end}-repmat(mean(new_x{end},2),1,size(new_x{end},2));
new_x{2} = proba_feature1';
% new_x{end} = NormalizeFea(new_x{end},0);
% new_x{end} = new_x{end}-repmat(mean(new_x{end},2),1,size(new_x{end},2));
end

