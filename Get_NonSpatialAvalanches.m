function [Size,Duration] = Get_NonSpatialAvalanches(Pop,th)
% Neuronal avalanches detected using a spatially unconstraint definition 
% (e.g., Beggs and Plenz). In this definition, the sum activity across all
% N recorded neurons is calculated at each time bin. Then, an avalanche is
% defined as a period when this summed population activity exceeds 
% a threshold th (nb. of neurons).
%
% Inputs:
% - Pop : sum activity
% - th : threshold (nb. of neurons)
%
% Outputs:
% - Size : size of avalanches
% - Duration : duration of avalanches
%
% Ponce-Alvarez A. 22/01/2024
%--------------------------------------------------------------------------
L = length(Pop);

Size = [];
Duration = [];

% Avalanches:
s = 0;
d = 0;
for i=2:L
    if Pop(i)>th
        s = s + Pop(i);
        d = d + 1;
    elseif Pop(i)<th && Pop(i-1)>th
        Size = [Size s];
        Duration = [Duration d];
        s = 0;
        d = 0;
    end
end
Size = Size(Duration>0);
Duration = Duration(Duration>0);
