function [P,K,L] = SL_Regulador_Robusto_alternativo(N,parms,mi,alfa)

%% ---------------------------------------------------

% Check input values

  if mi <= 0
     error('>> Valor de mi deve ser positivo.')
  end
 
   if (parms ~= 0) && (parms ~= 1)
      error('>> parms = {0,1}')
  end 
 
%% ---------------------------------------------------

% Global variables
  global m n q F G Q R H Ef Eg P0
  
%% ---------------------------------------------------

% Check weighting matrices  

  if min( eig( P0 ) ) < 0 
     error('P0 matrix is not positive semidefinite.') 
  end
  
  if min( eig( Q ) ) < 0 
     error('Q matrix is not positive semidefinite.') 
  end
  
  if min( eig( R ) ) <= 0 
     error('R matrix is not positive definite.') 
  end

%% --------------------------------------------------- 

% Auxiliary variables
  In = eye(n);
  Im = eye(m);
  Iq = eye(q);
  
%% ---------------------------------------------------- 

% Auxiliary matrices

  P(:,:,1) = P0;
  
  Qbar = [ zeros(n) In; In -Q ];
  Rbar = [ zeros(m) Im; Im -R ];
  
  A_cal = [ eye(n)      zeros(n,m); ...
            zeros(n)    zeros(n,m); ...
            zeros(m,n)  Im        ; ...
            zeros(m,n)  zeros(m)  ; ...
            zeros(n)    zeros(n,m); ...
            zeros(n)    zeros(n,m); ...
            In         -G         ; ...
            zeros(q,n) -Eg        ];
          
  % % posto coluna pleno
  % postoAc = [ size(A_cal,2) rank(A_cal) ]
  % pause       
  %%%--------------------------------------------------
  
  Z = [ zeros(5*n+2*m+q,n) ; In ; zeros(m,n) ];
    
  V = [ zeros(6*n+2*m+q,m) ; Im ];
    
  U = [ zeros(2*n,n) ; zeros(2*m,n) ; -In ; zeros(n); F ; Ef ; zeros(n+m,n) ];
    
%% ----------------------------------------------------

% Finding lambda
  lambda = (1 + alfa) * norm( mi * H' * H );
          
%%%=====================================================================%%%
%%%               CONDITION FOR THE EXISTENCE OF THE SOLUTION           %%%
%%%=====================================================================%%%
        
%%%---------------------------------------------------
%%% Condições sobre o posto   
%%%---------------------------------------------------
      
% Posto linha pleno
  Aux = [ In -G ; zeros(q,n) -Eg ];
   if ( rank( Aux ) < size( Aux,1 ) ) && parms ==0
      error(' No solution for PARMS = 0') 
  end
     
%%%---------------------------------------------------
%%% Conditions on Phi(mi,lambda)  
%%%---------------------------------------------------
    
phi = mi^(-1) * In - lambda^(-1) * H * H';
%-----------------------------------------
% Matrix phi(mi,lambda) >= 0 
   if min( eig( phi ) ) < 0
      error('ErrorTests:convertTest', ...
      '>> A matriz Phi não é (semi)definida positiva \n Modifique as variáveis parms e/ou mi.')
   end

   
%%%=====================================================================%%%
%%%                         Loop                                        %%%
%%%=====================================================================%%%

% Computing Sigma
  MSigma = parms * blkdiag( phi,lambda^(-1)*Iq );

%%%----------------------------------------------------

for k = 1:N-1
    
    Pbar = [ zeros(n) In ; In -P(:,:,k) ];
    
    P_cal = blkdiag( Pbar,Rbar,Qbar );
                  
    Aux_W = blkdiag( P_cal,MSigma );
            
    W = [ Aux_W  A_cal ; A_cal' zeros(n+m) ];
    
    %%% Check number and condiiton for matrix W 
    %%%---------------------------------------------------- 
    if rcond(W) <= eps
       error('ErrorTests:convertTest', ...
       'W is ill conditioned \n Chance parms and/or mi.')
    end   
    
    %%%=================================================================%%%
    %%%              [L; K; P] = [Z V U ]^T * inv(W) * U                %%%
    %%%=================================================================%%%
    
    Sol = [ Z V U ]' * inv(W) * U;
    
    AUX = mat2cell( Sol,[n m n],[n] );
    
      L(:,:,k) = AUX{1};
    
      K(:,:,k) = AUX{2};
    
    P(:,:,k+1) = AUX{3};
    
    %%%---------------------------------------------------    
    % Check if Riccati isn't positive semidefinite
      if min( eig( P(:,:,k+1) ) ) < 0
         error('>> A Riccati não é (semi)definida positiva.') 
      end
    %%%---------------------------------------------------
 end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%