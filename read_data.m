function [xs] = read_data(src, pt)
root_folder = '';

testdata=cell(1,9);

testdata{1} = importdata(strcat(src,pt,'_DP.txt'));
testdata{1} = NormalizeFea(testdata{1},1);
testdata{2} = importdata(strcat(src, pt ,'_DR.txt'));
testdata{2} = NormalizeFea(testdata{2},1);
testdata{3} = importdata(strcat(src , pt , '_DT.txt'));
testdata{3} = NormalizeFea(testdata{3},1);
testdata{4} = importdata(strcat(src , pt , '_PseAAC.txt'));
testdata{4} = NormalizeFea(testdata{4},1);
testdata{5} = importdata(strcat(src , pt , '_Top-n-gram.txt'));
testdata{5} = NormalizeFea(testdata{5},1);
testdata{6} = importdata(strcat(src , pt , '_Kmer.txt'));
testdata{6} = NormalizeFea(testdata{6},1);
testdata{7} = importdata(strcat(src , pt , '_PPCT.txt'));
testdata{7} = NormalizeFea(testdata{7},1);
testdata{8} = importdata(strcat(src , pt , '_PSFM_DBT.txt'));
testdata{8} = NormalizeFea(testdata{8},1);
testdata{9} = importdata(strcat(src , pt , '_PSSM_DT.txt'));
testdata{9} = NormalizeFea(testdata{9},1);
testdata{10} = importdata(strcat(src , pt , '_Bit20NTCT(2)40.txt'));
testdata{10} = NormalizeFea(testdata{10},1);
testdata{11} = importdata(strcat(src , pt , '_CTD3.txt'));
testdata{11} = NormalizeFea(testdata{11},1);
testdata{12} = importdata(strcat(src , pt , '_GAP5.txt'));
testdata{12} = NormalizeFea(testdata{12},1);
testdata{13} = importdata(strcat(src , pt , '_Ngram(N=1)1.txt'));
testdata{13} = NormalizeFea(testdata{13},1);
testdata{14} = importdata(strcat(src , pt , '_OVNT(5)84.txt'));
testdata{14} = NormalizeFea(testdata{14},1);

load([root_folder, 'matlab_models/', pt ,'/RF.mat'], 'RFModel')
load([root_folder, 'matlab_models/', pt ,'/M.mat'],'M')
load([root_folder, 'matlab_models/', pt ,'/mean_train.mat'],'mean_train')

ttdata = testdata;

for iv = 1:length(ttdata)
%    Center
    ttdata{iv} = ttdata{iv}'-repmat(mean_train{iv}, 1, size(ttdata{iv}',2));
end



new_xs = Proba_Class_feat_RF(RFModel, ttdata );

[new_xs,~]=ZscoreNormalize(new_xs,1,M);
new_xs = cellfun(@(x) x', new_xs, 'UniformOutput', false);

xs = [ttdata,new_xs];

nst = cellfun(@(x) x', xs, 'UniformOutput', false);

xs = nst;

end

