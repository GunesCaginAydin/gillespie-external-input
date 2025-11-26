function [W,Ampli,NE,NI,typ,xyz] = Get_Connectivity(conntype,varargin)

% Constructs connectivity
% Inputs:
% Connectivity type:
% conntype = 1 : all-to-all connectivity
% conntype = 2 : sparse
% conntype = 3 : spatially-embedded
%    If conntype = 3, provide decay parameter.
%
% Outputs:
% - W : connectivity matrix
% - Ampli : balanced amplification index
% - NE, NI : nb. of E and I neurons
% - type: cell-type array
%
%
% Adrian Ponce-Alvarez 11/07/2024
%--------------------------------------------------------------------------


switch conntype

    case 1 % All-to-all connectivity

    % Total number of neurons:
    N = 2000;
    n = N/2;

    w0 = 0.2;
    wE = 7;
    wI = wE-w0;

    WE = wE*ones(n);
    WI = wI*ones(n);

    W = [WE -WI;WE -WI];    

    NE = n;
    NI = n;
    typ = zeros(1,N);
    typ(1:NE) = 1;

    T = schur(W);
    [~,d] = eig(T);
    d = diag(d);
    Y = T*T';
    Ampli = 1 - sum( abs(d).^2 )/trace(Y);    
    xyz = [];

    case 2 % Sparse connectivity with balanced amplification

    % Total number of neurons:
    N = 2000;
    n = N/2;
        
    % sparsity:
    f1 = .1;
    f2 = .001;
    
    J = 10;
    
    q = rand(n,1)<f1;
    D1 = J*f1*diag(q);
    q = rand(n,1)<f1;
    D2 = J*f2*diag(q);

    % Random Orthogonal Matrix:
    [Q,R] = qr(randn(n));
    Q = Q*diag(sign(diag(R)));
    [M,R] = qr(randn(n));
    M = M*diag(sign(diag(R)));
    
    WE = 0.5*( Q\D1*Q + M\D2*M );
    WI = 0.5*( Q\D1*Q - M\D2*M );
    
    WE(WE<0)=0;
    WI(WI<0)=0;
    
    W = [WE -WI;WE -WI];
    
    NE = n;
    NI = n;
    typ = zeros(1,N);
    typ(1:NE) = 1;

    T = schur(W);
    [~,d] = eig(T);
    d = diag(d);
    Y = T*T';
    Ampli = 1 - sum( abs(d).^2 )/trace(Y);    
    xyz = [];

    case 3 % Spatially embeded connectivity

    file = 'PreparedData_fish_3';
    dir = 'C:\Toni\Adrian\Codi\Dades\';
    load([dir file],'Type','avg','xyz')
    
    XYZ = xyz(Type>-1,:);
    typ = Type(Type>-1);
    
    N = length(typ);
    NE = sum(typ==1);
    NI = sum(typ==0);
    
        if nargin > 1
        lambda = varargin{1}; % connectivity decay in microns
        else
        disp('Error: provide a value for lambda')
        end
    
        D = zeros(N);
        for i =1:N
            ri = XYZ(i,:);
            for j=1:N
                rj = XYZ(j,:);
                D(i,j) = sqrt( sum((ri-rj).^2) ); % Euclidean dist.
            end
        end

    WEE = rand(NE) < exp(-D(1:NE,1:NE)/lambda);
    WIE = rand(NI,NE) < exp(-D(NE+1:end,1:NE)/lambda);
    WEI = rand(NE,NI) < exp(-D(1:NE,NE+1:end)/lambda);
    WII = rand(NI) < exp(-D(NE+1:end,NE+1:end)/lambda);
        
    wE = 1.4; %10;
    wI = 1.6; %12;
    
    WEE_norm = WEE./repmat(sum(WEE,2)/wE,[1,size(WEE,2)]);
    WEI_norm = WEI./repmat(sum(WEI,2)/wI,[1,size(WEI,2)]);
    WIE_norm = WIE./repmat(sum(WIE,2)/wE,[1,size(WIE,2)]);
    WII_norm = WII./repmat(sum(WII,2)/wI,[1,size(WII,2)]);
    
    W = [WEE_norm -WEI_norm;WIE_norm -WII_norm];
    
    T = schur(W);
    [~,d] = eig(T);
    d = diag(d);
    Y = T*T';
    Ampli = 1 - sum( abs(d).^2 )/trace(Y);
    xyz = XYZ;

end

return