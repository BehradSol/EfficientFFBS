%% This file is distributed under BSD (simplified) license
%% Author: Behrad Soleimani <behrad@umd.edu>

function [m, Cov, P] = Filtering(y, Ao, Qo, Co, R, e, Bo, m0, V0)
    % This function implements the conventional Forward Filtering and 
    % Backward Smooting.
    
    % Inputs:
    % y  = observation matrix with dimension N_y * T
    % Ao = VAR coefficients; A is a p * 1 cell corresponding to each lag 
    %     such that each cell is an N_x * N_x matrix
    % Qo = source noise covariance matrix with dimension N_x * N_x
    % Co = linear mapping matrix between observations and sources with
    %     dimension N_y * N_x
    % R  = observation noise matrix with dimension N_y * N_y
    % e  = stimuli matrix with dimension N_e * T
    % Bo = stimuli coefficients matrix with dimension N_x * N_e
    % m0 = initial of mean vector
    % V0 = initial covariance matrix
    

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
    Ny = length(R);
    L = size(y);
    T = L(1,2);
    
    if nargin < 6
        m0 = zeros(NNx,1);
        V0 = zeros(NNx,NNx);
        Ne = 1;
        B = zeros(Nx , Ne);
        e = zeros(Ne , T);
    end
    
    if nargin < 8
        m0 = zeros(NNx,1);
        V0 = zeros(NNx,NNx);
    end
    
    
    L = size(Bo);
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
       
    B = zeros(NNx , Ne);
    B(1:Nx , :) = Bo;
    % ---------------------------------------------------------------------
    % Coventional Filtering 
   
    mF = zeros(NNx,T);
    SF = zeros(NNx, NNx, T);

    mLast = m0;
    SLast = V0;
    
    mTemp = zeros(NNx,1);
    STemp = zeros(NNx,NNx);

    for i = 1 : 1 : T
        
        mTemp(Nx+1:end , :) = mLast(1: NNx - Nx , :);
        
        mTemp(1:Nx , :) = A(1:Nx , :)*mLast + Bo * e(:,i);
        STemp = A*SLast*A' + Q;
        
        K = STemp(:,1:Nx)*Co'/( Co*STemp(1:Nx,1:Nx)*Co'+ R);
        
        mF(:,i) = mTemp + K*(y(:,i) - Co*mTemp(1:Nx , :));
        SF(:,:,i) = STemp - K * (Co*STemp(1:Nx,1:Nx)*Co'+ R) *K';

        mLast = mF(:,i);
        SLast = SF(:,:,i);     
    end
    
    mB = zeros(NNx,T);
    mB(:,end) = mF(:,end);
    
    SB = zeros(NNx,NNx,T);
    SB(:,:,end) = SF(:,:,end);
    
    S = zeros(NNx,NNx,T);

    for i = T-1 : -1 : 1  
        S(:,:,i) = SF(:,:,i)*A'/(A*SF(:,:,i)*A' + Q);

        mB(:,i) = mF(:,i) + S(:,:,i)*(mB(:,i+1) - A*mF(:,i) - B * e(:,i));
        
        SB(:,:,i) = SF(:,:,i) + S(:,:,i)*(SB(:,:,i+1)- ( A*SF(:,:,i)*A' + Q) )*S(:,:,i)';
    end

    m = mB((p-1)*Nx + 1 : p*Nx,:); %Mean vectors with dimension Nx * T
    Cov = SB((p-1)*Nx + 1 : p*Nx , (p-1)*Nx + 1 : p*Nx , :);  %Covariance matrices with dimension Nx * Nx * T
    
    % ---------------------------------------------------------------------
    % Generating cross-covariance terms
    P = [];
    
    % Before run this part, be careful! It might have some memory issues!
    % You may change this part for your own problem if it is large-scale!
    
%     Po = cell(T,T);
%     P  = cell(T,T);
% 
%     for i = 1 : T
%         Po{i,i} = SB(:,:,i);
%         P{i,i}  = SB(1:Nx , 1:Nx ,i);
%     end
%     
%     for i = T : -1 : 1
%         for j = T : -1 : i+1
%             Po{i,j} = Po{i+1,j}*S(:,:,i);
%             Po{j,i} = Po{i,j}';
%             
%             P{i,j} = Po{i,j}(1:Nx , 1:Nx);
%             P{j,i} = P{i,j}'; 
%         end
%     end
%     clear Po
end

