
clear all
addpath(genpath(pwd));%include all subfolders

% Construct connectivity:
%---------------------------------------------
% Connectivity type:
% conntype = 1 : all-to-all connectivity
% conntype = 2 : sparse
% conntype = 3 : spatially-embedded
conntype = 3;
lambda = 200;
[W,Ampli,NE,NI,typ,xyz] = Get_Connectivity(conntype,lambda);
N = NE + NI;

% Model parameters:
%---------------------------------------------
% transfer functions:
%response_fn = @(x) 1./(1+exp(-x));
response_fn = @(x) tanh(x).*(x>0);

% params:
beta_param = 1;
alpha_param = 0.1;

% time window:
t_min = 0;
t_max = 1800;
n_batch = 3;

% Inputs:
%I = 0.001*ones(N,1);
I = 0.001*ones(N,t_max*2);

% run simulation:
%---------------------------------------------
spike_times = [];
spike_ids   = [];
time = [];
Fds = [];
sh = 0;
for n = 1:n_batch

    % initial state:
    if n==1
    init_state = zeros(2,N);
    init_state(1,:) = rand(1,N)<.3;%TG: decisio de si son actives o passives a l'inici
    init_state(2,:) = ~init_state(1,:);%complementari
    else
    init_state = network_state;
    end
    
    [sp_times,sp_ids,network_state] = ...
        Gillespie_EImodel(W,response_fn,beta_param,alpha_param,I,t_min,t_max,init_state);

    shift = (n-1)*t_max;
    spike_times = [spike_times (sp_times + shift)];
    spike_ids   = [spike_ids sp_ids];

    % % Convolution -> fluorescence signal:
    %--------------------------------------
    tauR = .01; % rise time (in the order of ms)
    tauD = .5;  % decay time (in the order of 0.5-2 s)
    K = .6;  % between 0.5-3
    q = 5;  % between 1-10
    Fm = 10; % around 10
    
    resol_F = 1/15;
    [t,fds] = SpikesToFluoresence(sp_times,sp_ids,N,resol_F,tauR,tauD,K,q,Fm);
    time = [time (t+sh)];
    sh = time(end);
    Fds = [Fds;fds];
    clear t fds
end

Tmax = t_max+shift;

% E activity:
E_spike_times = spike_times(spike_ids<=NE);
E_spike_ids = spike_ids(spike_ids<=NE);

% I activity:
I_spike_times = spike_times(spike_ids>NE);
I_spike_ids = spike_ids(spike_ids>NE);

% Firing rates:
Rates = zeros(1,N);
for i=1:N
Rates(i) = sum(spike_times(spike_ids==i))/(Tmax);
end
MeanRate = mean(Rates);
firing_rates = Rates;


% sum spiking activity:
dT = 1; % (Shew et al.)
bins = 0:dT:Tmax;
nbins =length(bins)-1;
Pop = zeros(1,nbins);
PopE = zeros(1,nbins);
PopI = zeros(1,nbins);
for i = 1:nbins
    Pop(i) = sum( spike_times>bins(i) & spike_times<bins(i+1) )/N;
    PopE(i) = sum( E_spike_times>bins(i) & E_spike_times<bins(i+1) )/NE;
    PopI(i) = sum( I_spike_times>bins(i) & I_spike_times<bins(i+1) )/NI;
end
EIratio = PopE./(PopE + PopI);

% avalanches spikes:
%-------------------------------------------------------

X = Pop*N;
th = floor(0.005*N);
[Size,Duration] = Get_NonSpatialAvalanches(X,th);
Nav = length(Size);

% <S>(T) function:
w = logspace(0,1.1*log10(max(Duration)),10);
Ts = nan(1,length(w)-1);
S = nan(1,length(w)-1);
for k=1:length(w)-1  
    ii=find(Duration>=w(k) & Duration<w(k+1));
    dur = Duration(ii);
    s   = Size(ii);
    Ts(k) = mean(dur);
    S(k) = mean(s);
end


% Get exponents:
%-----------------------------
% Maximum Likelihood Estimation:
tau = plmle(Size,'xmin',10,'xmax',10*N);
Alpha = plmle(Duration,'xmin',min(Duration),'xmax',max(Duration));
% Least-squares for <S>(T):
cut = 1; %*resol; 
X = Ts(Ts>=cut);
Y = S(Ts>=cut);
logx = log(X);
logy = log(Y);
[~, a2, ~, Ea2, R_PL] = linear_fit([logx' logy'],min(logx),max(logx),1,0);
sigmaNuZ = 1/a2;


% Avalanches calcium events:
%----------------------------------

% Gaussian model:
%Raster_F = GetCalciumEvents(Fds,3);

%z-score:
SD = std(Fds);
F = (Fds - repmat(mean(Fds),[size(Fds,1),1]))./repmat(SD,[size(Fds,1),1]);
% Binary events:
Raster_F = F > 3;

% Firing rates:
Rates_F = sum(Raster_F)/(Tmax);
MeanRate_F = mean(Rates_F);

% Avalanches:
X = sum(Raster_F,2);
th = floor(0.005*N);
[Size_F,Duration_F] = Get_NonSpatialAvalanches(X,th);
Nav_F = length(Size_F);

% <S>(T) function:
w = logspace(0,1.1*log10(max(Duration_F)),15);
Ts_F = nan(1,length(w)-1);
S_F = nan(1,length(w)-1);
for k=1:length(w)-1  
    ii=find(Duration_F>=w(k) & Duration_F<w(k+1));
    dur = Duration_F(ii);
    s   = Size_F(ii);
    Ts_F(k) = mean(dur);
    S_F(k) = mean(s);
end

% Get exponents:
%-----------------------------
% Maximum Likelihood Estimation:
tau_F = plmle(Size_F,'xmin',10,'xmax',10*N);
Alpha_F = plmle(Duration_F,'xmin',min(Duration_F),'xmax',max(Duration_F));
% Least-squares for <S>(T):
cut = 1; %*resol; 
X = Ts_F(Ts_F>=cut);
Y = S_F(Ts_F>=cut);
logx = log(X);
logy = log(Y);
[~, a2, ~, Ea2, R_PL] = linear_fit([logx' logy'],min(logx),max(logx),1,0);
sigmaNuZ_F = 1/a2;

% Spatial definition (Ponce-Alvarez et al. 2018):
%----------------------------------------------------------------------
[Size_Space,Duration_Space] = avalanche_analysis_EI(Raster_F,xyz,typ,10,3);

% <S>(T) function:
w = logspace(0,1.1*log10(max(Duration_Space)),15);
Ts_Space = nan(1,length(w)-1);
S_Space = nan(1,length(w)-1);
for k=1:length(w)-1  
    ii=find(Duration_Space>=w(k) & Duration_Space<w(k+1));
    dur = Duration_Space(ii);
    s   = Size_Space(ii);
    Ts_Space(k) = mean(dur);
    S_Space(k) = mean(s);
end

% Get exponents:
%-----------------------------
% Maximum Likelihood Estimation:
tau_Space = plmle(Size_Space,'xmin',th,'xmax',10000);
Alpha_Space = plmle(Duration_Space,'xmin',5,'xmax',30/resol_F);
% Least-squares for <S>(T):
cut = 4; %*resol; 
X = Ts_Space(Ts_Space>=cut);
Y = S_Space(Ts_Space>=cut);
logx = log(X);
logy = log(Y);
[~, a2, ~, Ea2, R_PL] = linear_fit([logx' logy'],min(logx),max(logx),1,0);
sigmaNuZ_Space = 1/a2;



% Figures:
%-----------------------------

figure
xSize = 17; ySize = 12.5;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[50 50 xSize*50 ySize*50],'Color','w')

axes('Position',[.065 .8 .22 .18])
x = 0:1000;
Pr = exp(-x/lambda);
area(x,Pr,'FaceColor',[.3 .7 .88],'EdgeColor',[.2 .4 .6],'LineWidth',1)
xlabel('Distance \itr\rm [\mum]','fontsize',9)
ylabel('Connection prob.','fontsize',9)
text(.29,.6,'P(\itr\rm) = exp(-\itr\rm/\lambda)','units','normalized','fontsize',10)
box off
text(-.2,1.01,'A','fontsize',12,'units','normalized','fontweight','bold')

axes('Position',[.065 .4 .22 .3])
% Raster: spiking activty:
t_lim = 10*60; % 10 min
ii = find(spike_times<t_max);
spT = spike_times(ii);
spID = spike_ids(ii);
col = lines(2);
hold on
plot(spT(spID<=NE),spID(spID<=NE),'.','markersize',3,'color',col(1,:))
plot(spT(spID>NE),spID(spID>NE),'.','markersize',3,'color',col(2,:))
set(gca,'xlim',[t_lim 2*t_lim],'YLim',[0 N],'xtick',[],'xcolor','w','fontsize',9)
box off
xlabel('Time [s]','fontsize',9)
ylabel('neuron ID','fontsize',9)
text(.32,.98,sprintf(' \\lambda = %g \\mum \n \\phi = %2.2f',lambda,Ampli),'units','normalized','EdgeColor','k','BackgroundColor','w')
box off
text(-.2,1.02,'B','fontsize',12,'units','normalized','fontweight','bold')

axes('Position',[.065 .26 .22 .11])
% Population activity:
hold on
plot(bins(1:end-1),PopE,'color',col(1,:))
plot(bins(1:end-1),PopI,'color',col(2,:))
set(gca,'xlim',[t_lim 2*t_lim],'fontsize',9)
box off
xlabel('Time [s]','fontsize',9)
ylabel('Pop. activity','fontsize',9)

axes('Position',[.065 .08 .22 .11])
% Population activity fluorescence:
Eact = sum(Raster_F(:,1:NE),2)/NE;
Iact = sum(Raster_F(:,NE+1:end),2)/NI;
hold on
plot(time(time<t_lim),Eact(time<t_lim),'color',col(1,:))
plot(time(time<t_lim),Iact(time<t_lim),'color',col(2,:))
set(gca,'xlim',[t_lim 2*t_lim],'fontsize',9)
box off
xlabel('Time [s]','fontsize',9)
ylabel('Pop. activity','fontsize',9)


axes('Position',[.38 .745 .16 .195])
% avalanches spiking activity
% Size distribution:
hold on
Bins = logspace(log10(min(Size)),1.1*log10(max(Size)),12);
[x,n]=get_pdfbins(Size,Bins);
plot(x,n,'o-','color',lines(1),'markerfacecolor',lines(1),'markersize',3)    
set(gca,'xscale','log','yscale','log','fontsize',9)
xlabel('Avalanche size S','fontsize',9)
ylabel('Probability density','fontsize',9)

xlim = [.8*x(find(~isnan(x),1,'first')) 1.2*x(find(~isnan(x),1,'last'))];

xo = xlim;
y = xo.^(-tau);
plot(xo,y,'--','color','k','linewidth',1)
set(gca,'xlim',xlim)
text(.1,.26,'P(S) \sim S^{-\tau}','units','normalized','fontsize',9)
text(.1,.1,sprintf('\\tau = %2.2f',tau),'units','normalized','fontsize',8)

box off
text(-.27,1.05,'C','fontsize',12,'units','normalized','fontweight','bold')


axes('Position',[.60 .745 .16 .195])
% avalanches spiking activity
% Duration distribution:
hold on
Bins = logspace(log10(min(Duration)),log10(max(Duration)),9);
[x,n]=get_pdfbins(Duration,Bins);
plot(x,n,'o-','color',lines(1),'markerfacecolor',lines(1),'markersize',3)    
set(gca,'xscale','log','yscale','log','fontsize',9)
xlabel('Avalanche duration T','fontsize',9)
ylabel('Probability density','fontsize',9)


xlim = [.8*x(find(~isnan(x),1,'first')) 1.2*x(find(~isnan(x),1,'last'))];
xo = xlim;
y = xo.^(-Alpha);
plot(xo,y,'--','color','k','linewidth',1)
set(gca,'xlim',xlim)
text(.1,.26,'P(T) \sim T^{-\alpha}','units','normalized','fontsize',9)
text(.1,.1,sprintf('\\alpha = %2.2f',Alpha),'units','normalized','fontsize',8)

box off
text(-.27,1.05,'D','fontsize',12,'units','normalized','fontweight','bold')



axes('Position',[.82 .745 .16 .195])
% avalanches spiking activity
% <S>(T):
plot(Ts,S,'s-','color',lines(1),'markerfacecolor',lines(1),'markersize',3)
set(gca,'xscale','log','yscale','log')

xlabel('Avalanche duration T','fontsize',9)
ylabel('<S>(T)','fontsize',9)

text(.4,.26,'<S>(T) \sim T^{1/\sigma\nuz}','units','normalized','fontsize',9)
text(.4,.1,sprintf('\\sigma\\nuz = %2.2f',sigmaNuZ),'units','normalized','fontsize',8)

box off

text(-.27,1.05,'E','fontsize',12,'units','normalized','fontweight','bold')



axes('Position',[.38 .41 .16 .195])
% avalanches fluorescence
% Size distribution:
hold on
Bins = logspace(log10(min(Size_F)),1.1*log10(max(Size_F)),12);
[x,n]=get_pdfbins(Size_F,Bins);
plot(x,n,'o-','color',lines(1),'markerfacecolor',lines(1),'markersize',3)    
set(gca,'xscale','log','yscale','log','fontsize',9)
xlabel('Avalanche size S','fontsize',9)
ylabel('Probability density','fontsize',9)


xlim = [.8*x(find(~isnan(x),1,'first')) 1.2*x(find(~isnan(x),1,'last'))];

xo = xlim;
y = xo.^(-tau_F);
plot(xo,y,'--','color','k','linewidth',1)
set(gca,'xlim',xlim)
text(.1,.26,'P(S) \sim S^{-\tau}','units','normalized','fontsize',9)
text(.1,.1,sprintf('\\tau = %2.2f',tau_F),'units','normalized','fontsize',8)

box off
text(-.27,1.05,'F','fontsize',12,'units','normalized','fontweight','bold')


axes('Position',[.60 .41 .16 .195])
% avalanches fluorescence
% Duration distribution:
hold on
Bins = logspace(log10(min(Duration_F)),log10(max(Duration_F)),12);
[x,n]=get_pdfbins(Duration_F,Bins);
plot(x,n,'o-','color',lines(1),'markerfacecolor',lines(1),'markersize',3)    
set(gca,'xscale','log','yscale','log','fontsize',9)
xlabel('Avalanche duration T','fontsize',9)
ylabel('Probability density','fontsize',9)


xlim = [.8*x(find(~isnan(x),1,'first')) 1.2*x(find(~isnan(x),1,'last'))];
xo = xlim;
y = xo.^(-Alpha_F);
plot(xo,y,'--','color','k','linewidth',1)
set(gca,'xlim',xlim)
text(.1,.26,'P(T) \sim T^{-\alpha}','units','normalized','fontsize',9)
text(.1,.1,sprintf('\\alpha = %2.2f',Alpha_F),'units','normalized','fontsize',8)

box off
text(-.27,1.05,'G','fontsize',12,'units','normalized','fontweight','bold')


axes('Position',[.82 .41 .16 .195])
% avalanches fluorescence
% <S>(T):
plot(Ts_F,S_F,'s-','color',lines(1),'markerfacecolor',lines(1),'markersize',3)
set(gca,'xscale','log','yscale','log')

xlabel('Avalanche duration T','fontsize',9)
ylabel('<S>(T)','fontsize',9)

text(.4,.26,'<S>(T) \sim T^{1/\sigma\nuz}','units','normalized','fontsize',9)
text(.4,.1,sprintf('\\sigma\\nuz = %2.2f',sigmaNuZ_F),'units','normalized','fontsize',8)

box off
text(-.27,1.05,'H','fontsize',12,'units','normalized','fontweight','bold')



axes('Position',[.38 .07 .16 .195])
% avalanches fluorescence
% Spatial definition
% Size distribution:
hold on
Bins = logspace(log10(min(Size_Space)),1.1*log10(max(Size_Space)),10);
[x,n]=get_pdfbins(Size_Space,Bins);
plot(x,n,'o-','color',lines(1),'markerfacecolor',lines(1),'markersize',3)    
set(gca,'xscale','log','yscale','log','fontsize',9)
xlabel('Avalanche size S','fontsize',9)
ylabel('Probability density','fontsize',9)


xlim = [.8*x(find(~isnan(x),1,'first')) 1.2*x(find(~isnan(x),1,'last'))];

xo = xlim;
y = xo.^(-tau_Space);
plot(xo,y,'--','color','k','linewidth',1)
set(gca,'xlim',xlim)
text(.1,.26,'P(S) \sim S^{-\tau}','units','normalized','fontsize',9)
text(.1,.1,sprintf('\\tau = %2.2f',tau_Space),'units','normalized','fontsize',8)

box off
text(-.27,1.05,'I','fontsize',12,'units','normalized','fontweight','bold')


axes('Position',[.60 .07 .16 .195])
% avalanches fluorescence
% Spatial definition
% Duration distribution:
hold on
Bins = logspace(log10(min(Duration_Space)),log10(max(Duration_Space)),13);
[x,n]=get_pdfbins(Duration_Space,Bins);
plot(x,n,'o-','color',lines(1),'markerfacecolor',lines(1),'markersize',3)    
set(gca,'xscale','log','yscale','log','fontsize',9)
xlabel('Avalanche duration T','fontsize',9)
ylabel('Probability density','fontsize',9)


xlim = [.8*x(find(~isnan(x),1,'first')) 1.2*x(find(~isnan(x),1,'last'))];
xo = xlim;
y = xo.^(-Alpha_Space);
plot(xo,y,'--','color','k','linewidth',1)
set(gca,'xlim',xlim)
text(.1,.26,'P(T) \sim T^{-\alpha}','units','normalized','fontsize',9)
text(.1,.1,sprintf('\\alpha = %2.2f',Alpha_Space),'units','normalized','fontsize',8)

box off
text(-.27,1.05,'J','fontsize',12,'units','normalized','fontweight','bold')


axes('Position',[.82 .07 .16 .195])
% avalanches fluorescence
% Spatial definition
% <S>(T):
plot(Ts_Space,S_Space,'s-','color',lines(1),'markerfacecolor',lines(1),'markersize',3)
set(gca,'xscale','log','yscale','log','xlim',[3 50])

xlabel('Avalanche duration T','fontsize',9)
ylabel('<S>(T)','fontsize',9)

text(.4,.26,'<S>(T) \sim T^{1/\sigma\nuz}','units','normalized','fontsize',9)
text(.4,.1,sprintf('\\sigma\\nuz = %2.2f',sigmaNuZ_Space),'units','normalized','fontsize',8)

box off
text(-.27,1.05,'K','fontsize',12,'units','normalized','fontweight','bold')


annotation('textbox',[0.5 .999 0.4 .01],'string','Non-spatial avalanches: spiking activity','fontsize',12,'edgecolor','none','fontweight','bold')
annotation('textbox',[0.5 .669 0.4 .01],'string','Non-spatial avalanches: calcium events','fontsize',12,'edgecolor','none','fontweight','bold')
annotation('textbox',[0.52 .320 0.4 .01],'string','Spatial avalanches: calcium events','fontsize',12,'edgecolor','none','fontweight','bold')


dir3='C:\Users\Gunes\projects\gillespie-external-input\Grafics\';
figname = ['model_avalanches_lambda_new' num2str(lambda)];
exportgraphics(gcf,[dir3 figname '.tiff'],'Resolution',600)




