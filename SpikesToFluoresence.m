function [time,Fds] = SpikesToFluoresence(spike_times,spike_ids,N,resol_F,varargin)

% This function implements a phenomenological model that converts spike times 
% to synthetic fluorescence time series; as in Wei et al. (2020).
%
% Inputs:
% - spike_times : spike times (in s)
% - spike_ids : neuron index associated to each spike time
% - N : nb. of neurons
% - resol_F : temporal resolution of the resulting fluoresence signals (s)
% 
% Outputs:
% - time : time array
% - Fds : fluoresence signals
%
% Optional parameters:
% [time,Fds] = SpikesToFluoresence(spike_times,spike_ids,resol_F,...
%              tauR,tauD,K,q,Fm)
% - tauR : rise time (in the order of ms) | default : .01 s
% - tauD : decay time (in the order of 0.3-3 s) | default : .5 s
% - K : gain of the sigmoid function (between .5-3) | default : 0.6
% - q : half-activation param. (between 1-10) | default = 5;
% - Fm : maximum amplittude (around 10) | default : 10
%
% Reference:
% Wei Z, Lin B-J, Chen T-W, Daie K, Svoboda K, Druckmann S (2020) 
% A comparison of neuronal population dynamics measured with 
% calcium imaging and electrophysiology. PLoS Comput Biol 16(9): e1008198. 
% https://doi.org/10.1371/journal.pcbi.1008198
%
% A. Ponce-Alvarez 14/01/2024
%--------------------------------------------------------------------------


%N = max(spike_ids);
T = max(spike_times);

% spikes at times are converted to a latent variable, c(t), 
% by convolution with a double-exponential kernel:

% parameters:
if nargin > 4
tauR = varargin{1};
tauD = varargin{2};    
K    = varargin{3};
q    = varargin{4};
Fm   = varargin{5};
else % default values:
tauR = .01; % rise time (in the order of ms)
tauD = .5;  % decay time (in the order of few s)
K = .6;  % between .5-3
q = 5;  % between 1-10
Fm = 10; % around 10
end

dt = .01; % time bin
t = (dt:dt:T);
L = length(t);


c = zeros(N,L);

for n= 1:N
    tk = spike_times(spike_ids == n);
    for k = 1:length(tk)
    y =  (t>tk(k)).*exp( -(t-tk(k))/tauD );
    z =  (t>tk(k)).*( 1 - exp( -(t-tk(k))/tauR ) );
    x = y.*z;
    x(isnan(x)) = 0;
    c(n,:) = c(n,:) + x;
    end
end

% add noise:
sigma = 0.1*mean(c(:)); % Gaussian internal noise (10%)
c = c + sigma*randn(N,L);
% c(t) is truncated at zero if noise drove it to negative values:
c(c<0) = 0;

% Second, c(t) is converted to a synthetic fluorescence signal 
% through a sigmoidal function:

S = @(x) 1./( 1 + exp(-K*(x-q)) );

% fluorescence signal:
Y = Fm*feval(S,c');
sigma2 = 0.1*mean(Y(:)); % Gaussian external noise (10%)
F = Y + sigma2*randn(L,N);

% downsampling as data:
%resol_F = 1/15;
tds = 0:resol_F:T;
L = length(0:resol_F:T);
Fds = zeros(L-1,N);
for j = 1:L-1
   ii = ( t>=tds(j) & t<tds(j+1) ); 
   Fds(j,:) = mean(F(ii,:));  
end

time = tds(1:end-1);

return