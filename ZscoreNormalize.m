function [XX, M] = ZscoreNormalize(X,dim,M)
%z-score ��׼��
% X: cell of n �� m matrcis, n: num of samples, m: num of colomns
%dim : 1 :normalize by column; 2: normalize by row
%M : a cell of size 2:
%       M{1}: means vector computed from training data;
%       M{2}: stds vector computed from training data

N_views = length(X);
if isempty(M)
    M = {};
    XX = {};
    for i = 1 : N_views
        x = X{i};
        if dim == 2
            x = x';
        end
        m = mean(x);
        %����ÿһ�е�������׼��
        s = std(x);
        s(s == 0) = eps;
        %��ÿ��Ԫ��ִ��(x - m )/ s����
        r = (x - m) ./ s;
%         r((x - m) < 1e-6) = 0;
        if dim == 2
            r = r';
        end
        XX{i} = r;
        M{i,1} = m;
        M{i,2} = s;
    end
else
    for i = 1 : N_views
        x = X{i};
        if dim == 2
            x = x';
        end
        m = M{i,1};
        s = M{i,2};
        r = (x - m) ./ s;
        if dim == 2
            r = r';
        end
        XX{i} = r;
    end
end



