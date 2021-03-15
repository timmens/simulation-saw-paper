function [regime, beta]=getregime(tht,threshold)
if nargin<2 || isempty(threshold)
    threshold=0;
end
[T,p]=size(tht);
ntht = norms(tht,2,2);
ind = (2:T)';
regime = [1;ind(ntht(2:T)>threshold);T+1];

if nargout>1
    if threshold > 0
        tht(ntht <= threshold)=0;
    end
    beta = cumsum(tht);
end

