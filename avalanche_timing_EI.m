function [Shapes,avg_shapes,Shapes_EI,avg_shapes_EI,AvD,numAv] = ...
    avalanche_timing_EI(Avalanche,Av_cluster_sizes,Av_cluster_E,Av_cluster_I)

% This function computes the initiation and termination time of avalanches, 
% together with the shape of the avalanches.
% 
% Inputs: 
% - Avalanche (returned by avalanche_analysis.m): this tensor contains the 
% information about the propagating clusters, i.e.:
%     - rows: index of the cluster
%     - columns: time frames
% - AvSize : sizes of avalanches
%
% Outputs:
% - Shapes : shapes of the avalanches grouped by duration
% - avg_shapes : average shapes of the avalanches grouped by duration
% - Shapes_EI and avg_shapes_EI : same as Shapes and avg_shapes but
% separating E and I activity
% - AvD: 6-by-numAv matrix, where numAv is the total number of avalanches
%        - AvD(1,:): initiation time of each avalanche
%        - AvD(2,:): termination time of each avalanche
%        - AvD(3,:): size of each avalanche
%        - AvD(4,:): size of each avalanche but only E activity
%        - AvD(5,:): size of each avalanche but only I activity
%        - AvD(6,:): average E/I ratio within the avalanches

% Reference: Ponce-Alvarez et al. (2018) Neuron
%
% Adrián Ponce-Alvarez
% 07/03/2023
%--------------------------------------------------------------------------

 [Nrows,T] = size(Avalanche);
 % Total number of avalanches:
 av_count = 0;
 for row = 1:Nrows       
     Av = Avalanche(row,:);
     for t=1:T-1
         if Av(t)==1 && Av(t+1)==0
         av_count = av_count + 1;
         end
     end        
 end
    
    
    % initiation, termination, and size of avalanches:
    %----------------------------------------------------------------------
    numAv = av_count; % total number of avalanches
    Time = cell(1,numAv); % time since avalanche start  
    AvD  = zeros(3,numAv);
    Npar = cell(1,numAv);
    Epar = cell(1,numAv);
    Ipar = cell(1,numAv);
    
    av_count = 0;
    
    for row = 1:Nrows
       
        Av = Avalanche(row,:);
        count_t  = 0;
        t_run = [];
        n_run = [];
        E_run = [];
        I_run = [];
        tt    = [];
        
        for t=1:T-1
            if Av(t)==1 && Av(t+1)==1 % if an avalanche starts or continues... accumulate info
                count_t  = count_t + 1;
                t_run = [t_run count_t]; % time since avalanche start
                n_run = [n_run Av_cluster_sizes(row,t)]; % avalanche profile
                E_run = [E_run Av_cluster_E(row,t)]; % avalanche E profile
                I_run = [I_run Av_cluster_I(row,t)]; % avalanche I profile
                tt    = [tt t];
            elseif Av(t)==1 && Av(t+1)==0 % if an avalanche terminates... store variables
                count_t  = count_t + 1;
                t_run = [t_run count_t];
                n_run = [n_run Av_cluster_sizes(row,t)]; 
                E_run = [E_run Av_cluster_E(row,t)]; % avalanche E profile
                I_run = [I_run Av_cluster_I(row,t)]; % avalanche I profile                
                tt    = [tt t];
                av_count = av_count + 1;
                Time{av_count} = t_run; % time since avalanche start
                Npar{av_count} = n_run;
                Epar{av_count} = E_run;
                Ipar{av_count} = I_run;
                AvD(1,av_count) = tt(1);      % time of avalanche intiation
                AvD(2,av_count) = tt(end);    % time of avalanche termination
                AvD(3,av_count) = sum(n_run); % avalanche size
                AvD(4,av_count) = sum(E_run); % avalanche E size
                AvD(5,av_count) = sum(I_run); % avalanche I size
                AvD(6,av_count) = nanmean(E_run./(E_run + I_run)); % average E/I ratio
                
                % reset for next avalanche
                count_t  =0; 
                t_run = [];
                n_run = [];
                E_run = [];
                I_run = [];
                tt    = [];
            end
        end
        
    end
            
    % Avalanche shapes:
    Duration = cellfun(@numel,Time); % durations of avalanches
    Shapes = cell(max(Duration),1); % avalanche shapes
    avg_shapes = cell(max(Duration),1); % avalanche shapes

    Shapes_EI = cell(max(Duration),2); % avalanche shapes for E neurons
    avg_shapes_EI = cell(max(Duration),2); % avalanche shapes for I neurons
    
    for w = 1:max(Duration)
       
        indx = find( Duration==w );       
        
        if ~isempty(indx)      
            y  = zeros(length(indx),w);
            yE = zeros(length(indx),w);
            yI = zeros(length(indx),w);
            for i=1:length(indx)               
                n = Npar{indx(i)};
                y(i,:) = n;         
                nE = Epar{indx(i)};
                yE(i,:) = nE;               
                nI = Ipar{indx(i)};
                yI(i,:) = nI;                               
            end    
            Shapes{w} = [Shapes{w};y];
            Shapes_EI{w,1} = [Shapes_EI{w,1};yE];
            Shapes_EI{w,2} = [Shapes_EI{w,2};yI];
        end
        
        avg_shapes{w} = mean(Shapes{w});
        avg_shapes_EI{w,1} = mean(Shapes_EI{w,1});
        avg_shapes_EI{w,2} = mean(Shapes_EI{w,2});
    end
    

    return
    