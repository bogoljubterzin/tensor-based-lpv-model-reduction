function POD = TensorSVD(W,GRID,degrees,type);
%%%
%%% SYNTAX
%%%     POD = TensorSVD(W,GRID,degrees,type);
%%%
%%% DESCRIPTION
%%%    For a tensor W with rows indexed as (multiple) spatial, ...
%%%    
%%% INPUTS
%%%    W: tensor of order N
%%%
%%%    GRID: structured object that defines spatial and temporal grid.
%%%       See gridding.
%%%
%%%    DEGREES: row vector with N entries of approximation degrees. With 
%%%    type='hosvd'
%%%
%%%    TYPE: string defining the type of POD basis:
%%%         'hosvd'   for higher order tensor SVD method
%%%          
%%%  
%%%  OUTPUTS
%%%      POD:         structured data object, currently with fields
%%%      POD.Phi:     basis functions grouped in q block-rows of degrees(i)
%%%                   columns
%%%      POD.Phidot:  derivative of basis functions
%%%      POD.Phiddot: second derivative of basis functions
%%%      POD.svalues: singular values
%%%      POD.energy:  percentage of total data energy captured in basis
%%%      POD.degrees: vector of dimensions of basis
%%%
%%%  COPYRIGHTS AND OWNERSHIP
%%%    Siep Weiland
%%%    Created: May 28, 2008
%%%    Last changes: November, 2013


    W_tensor = tensor(W);

    %   Determine tensor order first
    N     = size(size(W_tensor),2);
    
% Determine number of different physical variables
    dimension  =  size(W_tensor);
        
    Phi     = zeros(max(GRID.dimension)+1,sum(degrees)); 
    Phidot  = zeros(max(GRID.dimension)+1,sum(degrees)); 
    Phiddot = zeros(max(GRID.dimension)+1,sum(degrees)); 
    svalues = zeros(sum(degrees),1);
    energy  = zeros(1,N);
    
   
    switch type
              
       case 'hosvd'
           disp('    ... doing POD basis computation with higher order SVD method ...  ');
       
           % Set up tensor first
           for i=1:N
               % make i-th unfolding
               Pi      =  tenmat(W_tensor, i, 'bc');  % lay-out tensor on ith axis as wide matrix
               szPi    =  size(Pi);                   % define its matrix size
               P = double(Pi*Pi')./szPi(2);           % define covariance matrix
                          
               %Compute singular vectors
               [phiPs,lambdas] = svd(P);              % ordinary SVD on covariance
               Psvs = sqrt(diag(lambdas));            % singular values w.r.t. axis i.
           
               %scale basis to reflect continous time inner product
               Gi     =  diag(GRID.delta(i) * ones(1,GRID.dimension(i)+1));  % simplest: uniform distribution
               Ri     =  chol(Gi);
               phiPs = Ri\phiPs;                      % scaled on sqrt of gramian
               Penergy = sum(diag(lambdas));

               % Define relative energy in basis functions
               % sum sv^2(k) /energy <= tol_energy
    
               %Take dominant basis functions
               phiPs = phiPs(:,1:degrees(i));
               
               % Approximate derivatives of basis functions
               phiPsdot=[]; 
               for k=1:degrees(i),
                  phiPsdot=[phiPsdot gradient(phiPs(:,k),GRID.delta(i))];
               end;  

               phiPsddot=[];    
               for k=1:degrees(i),
                   phiPsddot=[phiPsddot gradient(phiPsdot(:,k),GRID.delta(i))];
               end;  
       
              % store degree(i) basis elements together with their
              % derivatives
               ri_index        = 1:(GRID.dimension(i) + 1);
               ci_index        = (sum(degrees(1:(i-1))) + 1) : sum(degrees(1:i));
               Phi(ri_index,ci_index)     =  phiPs;
               Phidot(ri_index,ci_index)  = phiPsdot;
               Phiddot(ri_index,ci_index) = phiPsddot;
               energy(i)  = Penergy;
               %svalues(ci_index,1) = Psvs;
       
               % clear variables
               Pi=[]; P=[]; Gi=[]; Ri=[]; Penergy=[]; Psvs=[];
               phiPs=[]; phiPsdot=[]; phiPsddot=[];
           end;
           
       
   otherwise 
       disp('    >>> Unknown method for POD basis computation <<<');
       disp('    >>> I believe I coded COSVD method elsewhere <<<');
       
end;         

   
   
% OUTPUTS
%--------------------------------------------------------------------------
POD.degrees = degrees;
POD.Phi     = Phi;
POD.Phidot  = Phidot;
POD.Phiddot = Phiddot;
POD.energy  = energy;
POD.svalues = svalues;
