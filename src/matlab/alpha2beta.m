function beta = alpha2beta(alpha,regime)
[mbar,p]=size(alpha);
T = regime(mbar+1)-1;
beta = zeros(T,p);
for i=1:mbar
    beta(regime(i):regime(i+1)-1,:)=repmat(alpha(i,:),regime(i+1)-regime(i),1);
end
