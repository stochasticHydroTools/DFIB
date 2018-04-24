% ----------------------------------------------------------------------- %
% Author: Yuanxun Bill Bao
% Date: Oct 2013
% Description: velocity interpolation for the new divergence-free IB method
% same as spreadMAC2Dnew.m, except phiweights replaced by mex functions to 
% speed up computation
%
% Inputs:
% a - (vector) potential function, but a scalar field in 2D
% X - Marker locations [x,y].
% N - Domain size [Nx,Ny].
% h - Cartesian grid spacing.
% kID   - kernels 'stnd3pt', 'stnd4pt', etc.
% kID_d - derivative of kernels 'stnd3pt_d', 'stnd4pt_d' , etc.
% K - Parameter for flexible kernels.
% ----------------------------------------------------------------------- %

function U = interpMAC2Dnew(a,X,N,h,kID,kID_d,K)

Nx = N(1); Ny = N(2);
Nb = size(X,1);  
U  = zeros(Nb,2);

for k = 1:Nb
    s = [X(k,1)/h, X(k,2)/h];
    p = floor(s);
    r = s-p;
    j1 = mod((p(1)-2):(p(1)+3),Nx)+1;   
    j2 = mod((p(2)-2):(p(2)+3),Ny)+1;
    
   
    
    % 1st component
    %w2 = phiweights(r(1),kID,K)*phiweights(r(2),kID_d,K).';
    w2 = KernelGrid2D(r(1), kID  , r(2), kID_d, K);
    U(k,1)= -sum(sum(w2.*a(j1,j2)))/h; % grad delta requires dividing by h
    
    % 2nd component
    %w1 = phiweights(r(1),kID_d,K)*phiweights(r(2),kID,K).';
    w1 = KernelGrid2D(r(1), kID_d, r(2), kID  , K);
    U(k,2)=  sum(sum(w1.*a(j1,j2)))/h;
end