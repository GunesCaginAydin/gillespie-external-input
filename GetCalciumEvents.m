function [Raster_F] = GetCalciumEvents(Fds,L)
% Gaussian model:
% We calculate the ROI's baseline noise by fitting a
% gaussian to the distribution of negative values of deltaFoF
% Then detect the events larger than a multiple L of the Gaussian SD
% 
% Uses: mygaussfit.m
%--------------------------------------------------------------------------

[numFrames,N] = size(Fds);
sd_n = zeros(1,N);
mu_n = zeros(1,N);
% Gaussian fit:
for n= 1:N
    X = Fds(:,n);
    [smoothDist,x] = ksdensity(X);
    [~,indPeak]=max(smoothDist);
    xFit = x(1:indPeak);
    dataToFit=smoothDist(1:indPeak)/numFrames;
    [s,m]=mygaussfit(xFit',dataToFit);
    sd_n(n) = s;
    mu_n(n) = m;
end
%z-score:
F = (Fds - repmat(mu_n,[numFrames,1]))./repmat(sd_n,[numFrames,1]);
% Binary events:
Raster_F = F > L;

return