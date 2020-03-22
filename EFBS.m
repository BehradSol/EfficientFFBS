%% This file is distributed under BSD (simplified) license
%% Author: Behrad Soleimani <behrad@umd.edu>

function [m, Cov, P] = EFBS(y, Ao, Qo, C, R, e, B, m0)
    % The function implements the efficient Forward Filtering and Backward
    % Smooting.
    
    % Inputs:
    % y  = observation matrix with dimension N_y * T
    % Ao = VAR coefficients; A is a p * 1 cell corresponding to each lag 
    %     such that each cell is an N_x * N_x matrix
    % Qo = source noise covariance matrix with dimension N_x * N_x
    % C  = linear mapping matrix between observations and sources with
    %     dimension N_y * N_x
    % R  = observation noise matrix with dimension N_y * N_y
    % e  = stimuli matrix with dimension N_e * T
    % B  = stimuli coefficients matrix with dimension N_x * N_e
    % m0 = initial condition of mean vector

    % Outputs:
    % m   = estimated mean matrix with dimension N_x * T
    % Cov = estimated steady-state covariance matrix with dimension N_x *
    %       N_x
    % P   = cross-covariance terms; P is a T * T cell corresponding to each 
    %       time-pari such that each cell is an N_x * N_x matrix  
    
    % ---------------------------------------------------------------------
    % This part generates the augmneted model
    smallvalue = 1e-15;
    p = length(Ao);
    Nx = length(Ao{1});
    NNx = p*Nx;
    L = size(y);
    T = L(1,2);
        
    if nargin < 6
        m0 = zeros(NNx,1);
        Ne = 1;
        B = zeros(Nx , Ne);
        e = zeros(Ne , T);
    end
    
    if nargin < 8
        m0 = zeros(NNx,1);
    end
        
    L = size(B);
    Ne = L(1,2);
    
    A = [];
    temp = zeros((p-1)*Nx,p*Nx);
 
    for j = 1 : p
        A = [A, Ao{j}];
        for i = 1 : p-1
            if (i==j)
                temp((i-1)*Nx+1:(i)*Nx, (j-1)*Nx+1:(j)*Nx) = eye(Nx);
            end
        end
    end
    A = [A;temp];
    
   
    Q = Qo;
    for i = 1 : p-1
        Q = blkdiag(Q, smallvalue*eye(Nx));
    end
    

    % ---------------------------------------------------------------------
    % Efficient Forward Filtering and Backward Smoothing (EFBS)
   
    Sig = idare(A',A'*C',Q,C*Q*C' + R,Q*C',[]);   %Sigma_{k|k}

    if (isempty(Sig))
        disp('Augmented A matrix is not stable!')
    end
    
    SigP = A*Sig*A' + Q;    %Sigma_{k+1|k}

    K = SigP(:,1:Nx)*C'/( C*SigP(1:Nx,1:Nx)*C'+ R);   %Kalman gain
    S = Sig*A'/SigP;    %Smoothing gain
    
 
    mF = zeros(NNx,T);

    mLast = m0;
  
    mTemp = zeros(NNx,1);
    
    for i = 1 : 1 : T
        
        mTemp(Nx+1:end , :) = mLast(1: NNx - Nx , :);
        mTemp(1:Nx , :) = A(1:Nx , :)*mLast + B * e(:,i);
          
        mF(:,i) = mTemp + K*(y(:,i) - C*mTemp(1:Nx , :));

        mLast = mF(:,i);
             
    end
    
    
    mB = zeros(NNx,T);
    mB(:,end) = mF(:,end);
    
    Cov = zeros(p*Nx , p*Nx , T);
    Cov(: , : , T) = Sig;
    for i = T-1 : -1 : 1         
        mB(:,i) = mF(:,i) + S*(mB(:,i+1) - A*mF(:,i) - B * e(:,i));
        
        Cov(: , : , i) = Sig + S*(Cov(: , : , i+1) - SigP)*S';
    end
 
       
    m = mB( (p-1)*Nx + 1 : p*Nx,:); %Mean vectors with dimension Nx * T
    
    
    % ---------------------------------------------------------------------
    % Generating cross-covariance terms
    P = [];
    
    % Before run this part, be careful! It might have some memory issues!
    % You may change this part for your own problem if it is large-scale!
    
%     Po = cell(T,T);
%     P  = cell(T,T);
% 
%     for i = 1 : T
%         Po{i,i} = Cov(:,:,i);
%     end
%     
%     for i = T : -1 : 1
%         for j = T : -1 : i+1
%             Po{i,j} = Po{i+1,j}*S;
%             Po{j,i} = Po{i,j}';
%             
%             P{i,j} = Po{i,j}(1:Nx , 1:Nx);
%             P{j,i} = P{i,j}'; 
%         end
%     end
%     clear Po
    
end

