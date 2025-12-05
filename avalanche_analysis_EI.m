function [AvSize,Duration,AvSize_E,AvSize_I,Avalanche,Av_cluster_sizes,Av_Center,Av_cluster_E,Av_cluster_I] = ...
         avalanche_analysis_EI(raster,Centers,EI_labels,rpix,Dim)

% This code constructs the activity tensor,
% calculates the connected components at each time frame,
% detects the propagating clusters (neuronal avalanches), 
% and calculates the durations and sizes of avalanches.
% 
% Inputs: 
% - raster : binary activity N-by-T matrix, where N is the number of
% neurons and T is the number of time frames
% - Centers : location of neurons, NxDim matrix, where Dim is the dimension
% (2D or 3D)
% - EI_labels : identity of the neurons (1: excitatory, 0: inhibitory).
% - rpix : radius defining nearby neurons (in pixels)
% - Dim : dimension of the activity (3D or 2D)
%
% Outputs:
% - AvSize : sizes of avalanches
% - AvSize_E : sizes of avalanches but only E activity
% - AvSize_I : sizes of avalanches but only I activity
% - Duration: durations of avalanches
% - Avalanche: this binary matrix contains the information about the propagating
% clusters, i.e.:
%     - rows: index of the cluster
%     - columns: time frames
% - Av_cluster_sizes: size of propagating clusters:
%     - rows: index of the cluster
%     - columns: time frames     
% - Av_Center: this tensor tracks the center of mass of the propagating clusters,
%  i.e.:
%     - rows: index of the cluster
%     - columns: time frames
%     - 3rd dimension: (x,y,z) coordinates of the center of mass
%
% - Av_cluster_E and Av_cluster_E: same as Av_cluster_sizes but for E and I
% neurons

% Reference: Ponce-Alvarez et al. (2018) Neuron
%
% Adrián Ponce-Alvarez
% 14/03/2023
%--------------------------------------------------------------------------

[T,N] = size(raster);

% Construct tensor:
%-------------------------------------------------
% transform centers into indices of the 3D tensor:
Centers = floor(Centers);
x = Centers(:,1);
y = Centers(:,2);
z = Centers(:,3);

% Size of the 3D tensor:
Nx = max(x)+10;
Ny = max(y)+10;
Nz = max(z)+10;

% Map neurons to linear tensor's indices:
I = sub2ind([Nx,Ny,Nz],Centers(:,1),Centers(:,2),Centers(:,3));
LinIndx = zeros(N,2);
LinIndx(:,1) = I;
LinIndx(:,2) = 1:N;


% initialize matrices of size (nb of avalanches)x(time frames)
% a priori we don't know the nb of avalanches, thus for initialization we
% choose 2000.
Avalanche = zeros(2000,T);
Av_cluster_sizes   = zeros(2000,T);
Av_cluster_E       = zeros(2000,T);
Av_cluster_I       = zeros(2000,T);
Av_Center = zeros(2000,T,Dim);


for t = 1:T

X = raster(t,:); % activity at time t
Nactive = sum(X); % nb of active cells at time t
IDactive = find(X); % IDs of active cells at time t

% Construct activity tensor:
Tensor = zeros(Nx,Ny,Nz); % neuron centers
Tensor_coarse = zeros(Nx,Ny,Nz); % coarse

if mod(t,1000)==0
   fprintf('frame %g over %g\n',t,T) 
end
        
        % Fill the activity tensor:
        for n=1:Nactive

            neuron = IDactive(n); % active neurons in the current frame
            % coordinates (Coord is in pixels):
            xm = Centers(neuron,1);
            ym = Centers(neuron,2);
            zm = Centers(neuron,3);
            
            % Tensor entries are 1 for active neurons, 0 for inactive ones:
            Tensor(xm,ym,zm) = 1;

            % fill all neighboring pixels < R
            x = xm-rpix:xm+rpix;
            y = ym-rpix:ym+rpix;
            z = zm-rpix:zm+rpix;
            nx = length(x);
            ny = length(y);
            nz = length(z);
            u = zeros(nx*ny*nz,3);
            c = 0;
            for i=1:nx
                for j=1:ny
                    for k=1:nz
                    c = c+1;
                    u(c,:) = [x(i) y(j) z(k)];
                    end
                end
            end

            % check if they don't cross borders:
            x = u(:,1);
            y = u(:,2);
            z = u(:,3);
            u( x<=0 | y<=0 | z<=0 | x>Nx | y>Ny | z>Nz,:) = [];
            a = u(:,1)-xm;
            b = u(:,2)-ym;
            c = u(:,3)-zm;
            ind = find( ( a.^2 + b.^2 + c.^2 ) <= rpix^2);         
            
            % fill coarse
            pixs = u(ind,:);
            for i = 1:size(pixs,1)
            Tensor_coarse(pixs(i,1),pixs(i,2),pixs(i,3)) = 1;
            end
            clear u
        end
    
    
        % Get connected components(clusters):
        %----------------------------------------------------------------------
        cc = bwconncomp(Tensor_coarse); % connected comp. in the coarse-grained tensor        
        clear Tensor_coarse
        Nob = cc.NumObjects; % nuber of clusters
        Lin = find(Tensor);
        clear Tensor
        Comps = cell(1,Nob);
        Csize = zeros(1,Nob);
        for j=1:Nob
            pixs = cc.PixelIdxList{j};
            inter = intersect(pixs,Lin);
            Comps{j} = inter; % tensor indices of neurons composing cluster j
            Csize(j) = length(inter); % cluster size
        end
    
        % Keep only clusters with >2 neurons
        Comps = Comps( Csize>2 );
        nc = length(Comps);
        
        % From indices to neurons:
        % The first column of Ct indicates wether each neuron participates
        % in a cluster. The second column indicates in which cluster the
        % neuron participates.
        Ct = zeros(N,2);    
        for j=1:nc
            
            nn = Comps{j};     
            ii = find( ismember(LinIndx(:,1),nn) );
            units = LinIndx(ii,2);
            neurons = unique(units);           
            Ct(neurons,1) = 1;  % active (1) neuron
            Ct(neurons,2) = j;  % ... that belongs to cluster j
            
        end
    
    
        % Check if the components overlap with those of the previous time
        % frame
        %----------------------------------------------------------------------
        if t>1

            % Components up to now:
            occupied = unique( Ctprev(:,2) );
            occupied = occupied(occupied>0);
            non_occupied = setxor(1:2000,occupied); % non-occupied clusters
            counter = 0;
            
            temp = Ct(:,2);

            for j=1:nc    

                neuron_in_cluster = (temp==j) ; % neurons participating in cluster j at current time t    
                indx = neuron_in_cluster.*Ctprev(:,1); % neurons active at time t and t-1
                int = sum(indx)>0; 

                if int==1 % at least one neuron of this comp. was previously active 
                    prevclus = Ctprev(indx==1,2);                    
                    up = unique(prevclus);
                    if length(up)>1 % it can be that neurons participate in 2 clusters, if this is the case, the largest avalanche absorbs the small ones
                    si = histc( prevclus, up );
                    [~,index]=max(si);
                    Ct(neuron_in_cluster,2) = up(index);
                    else
                    Ct(neuron_in_cluster,2) = up;
                    end
                else % otherwise, new component 
                    counter = counter + 1;
                    Ct(neuron_in_cluster,2) = non_occupied(counter);
                end

            end
            
        end
        
        clus = unique( Ct(:,2) );
        clus = clus(clus>0); % current clusters
        
        % Avalanche : matrix of (cluster labels)x(time frames) size
        % with entries (k,t) equal to 1 if cluster k is active at time t
        Avalanche(clus,t) = 1;
        
        % Av_cluster_sizes : matrix of (cluster labels)x(time frames) size
        % with entries equal to the number of neurons active in cluster k
        % at time t
        %
        % Av_Center : matrix of (cluster labels)x(time frames) size
        % with entries equal to the center of mass of neurons active in cluster k
        % at time t        
        for j=1:length(clus)
           neurons = find( Ct(:,2) == clus(j) );  
           nb = sum( Ct(:,2) == clus(j) ); 
           Av_cluster_sizes(clus(j),t) = nb; % number of neurons in cluster
           Av_cluster_E(clus(j),t) = sum(EI_labels(neurons)==1); % number of E neurons in cluster
           Av_cluster_I(clus(j),t) = sum(EI_labels(neurons)==0); % number of I neurons in cluster
           cm = Centers(Ct(:,2) == clus(j),:); %centers of neurons in cluster j
           mcm = mean(cm); % center of mass of the avalnache
           Av_Center(clus(j),t,:) = mcm;
        end
        
        % update previous cluster information:
        Ctprev = Ct;


end

% clean matrices:
S = sum(Avalanche,2);
Avalanche(S==0,:)=[];
S = sum(Av_cluster_sizes,2);
Av_cluster_sizes(S==0,:)=[];
Av_cluster_E(S==0,:)=[];
Av_cluster_I(S==0,:)=[];
Av_Center(S==0,:,:)=[];

%--------------------------------------------------------------------------    
% matrices Avalanche, Av_cluster_sizes, and Av_Center summarize the avalanche activity
% Avalanche indicates which clusters were active at each time frame
% A_size indicates the nb of neurons participating in the clusters at each
% time
% Av_Center indicates the center of mass of the avalanche at each time
    
% Get avalanche sizes and durations
% avalanche size is defined as the total nb of activations during the
% avalanche (including neurons that activate more than one time during the
% avalanche)

    Life      = [];
    AvSize    = [];
    AvSize_E  = [];
    AvSize_I  = [];
    Nrows = size(Avalanche,1);
    for row = 1:Nrows
       
        Av = Avalanche(row,:);
        count_t  = 0;
        count_nb = 0;
        count_E  = 0;
        count_I  = 0;
        for t=1:T
            if Av(t)==1
                count_t  = count_t + 1;
                count_nb = count_nb + Av_cluster_sizes(row,t);
                count_E = count_E + Av_cluster_E(row,t);
                count_I = count_I + Av_cluster_I(row,t);
            else
                Life      = [Life count_t];
                AvSize    = [AvSize count_nb];
                AvSize_E  = [AvSize_E count_E];
                AvSize_I  = [AvSize_I count_I];
                count_t  =0; 
                count_nb =0; 
                count_E  =0; 
                count_I  =0; 
            end
        end
        
    end
    
    % Keep avalanches that last more than one time frame:
    Duration = Life(Life>1);
    AvSize   = AvSize(Life>1);
    AvSize_E = AvSize_E(Life>1);
    AvSize_I = AvSize_I(Life>1);
    
return    




