function [x,n,d,Ex,nraw]=get_pdfbins(X,bins)

S = length(X);
nbins = length(bins);
x = nan(1,nbins);
n = nan(1,nbins);
d = nan(1,nbins);
Ex = nan(1,nbins);
nraw = zeros(1,nbins);

bins = [bins bins(end)+bins(end)-bins(end-1)];

for w = 1:nbins
   
    y = X( X>=bins(w)-.0001 & X<bins(w+1)-.0001 );
    if ~isempty(y)
    x(w) = mean(y);
    n(w) = length(y)/S/(bins(w+1)-bins(w));
    Ex(w) = std(y); %/sqrt(length(y));
    nraw(w) = length(y);
    end
    d(w) = bins(w+1)-bins(w);
    
end



end