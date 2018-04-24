% ----------------------------------------------------------------------- %
% Author: Yuanxun Bill Bao
% Date: June 2014
% Description: velocity interpolation for the new divergence-free IB
% method, 3D problem, interpolation is done on the velocity potential
% function a(x,t), and a(x,t) is edge-centered.
%
% Inputs:
% a - velocity potential function
% X - Marker locations [x,y,z].
% N - Domain size [Nx,Ny,Nz].
% h - Cartesian grid spacing.
% kID   - kernels 'stnd3pt', 'stnd4pt', etc.
% kID_d - derivative of kernels 'stnd3pt_d', 'stnd4pt_d' , etc.
% K - Parameter for flexible kernels.
% ----------------------------------------------------------------------- %

function U = interpMAC3Dnew(a,X,N,h,kID,kID_d,K)

Nx = N(1); Ny = N(2); Nz = N(3);
Nb = size(X,1);
U  = zeros(Nb,3);

for k = 1:Nb
    
    % x1-edge
    s  = [X(k,1)/h-1/2, X(k,2)/h, X(k,3)/h];
    p  = floor(s);
    r1  = s-p;
    i1 = mod((p(1)-2):(p(1)+3),Nx)+1;
    i2 = mod((p(2)-2):(p(2)+3),Ny)+1;
    i3 = mod((p(3)-2):(p(3)+3),Nz)+1;
    
    % x2-edge
    s  = [X(k,1)/h, X(k,2)/h-1/2, X(k,3)/h];
    p  = floor(s);
    r2  = s-p;
    j1 = mod((p(1)-2):(p(1)+3),Nx)+1;
    j2 = mod((p(2)-2):(p(2)+3),Ny)+1;
    j3 = mod((p(3)-2):(p(3)+3),Nz)+1;
    
    % x3-edge
    s  = [X(k,1)/h, X(k,2)/h, X(k,3)/h-1/2];
    p  = floor(s);
    r3  = s-p;
    k1 = mod((p(1)-2):(p(1)+3),Nx)+1;
    k2 = mod((p(2)-2):(p(2)+3),Ny)+1;
    k3 = mod((p(3)-2):(p(3)+3),Nz)+1;
    
    
    % U1 = a2 * w3 - a3 * w2
    w3 = KernelGrid3D(r2(1),kID, r2(2),kID,   r2(3),kID_d, K);
    w2 = KernelGrid3D(r3(1),kID, r3(2),kID_d, r3(3),kID, K);
    U(k,1) = ( sum(sum(sum(w3.*a(j1,j2,j3,2)))) - sum(sum(sum(w2.*a(k1,k2,k3,3)))) ) / h;

    % U2 = a3 * w1 - a1 * w3
    w1 = KernelGrid3D(r3(1),kID_d, r3(2),kID,  r3(3),kID, K);
    w3 = KernelGrid3D(r1(1),kID,   r1(2),kID,  r1(3),kID_d, K);
    U(k,2) = ( sum(sum(sum(w1.*a(k1,k2,k3,3)))) - sum(sum(sum(w3.*a(i1,i2,i3,1)))) ) / h;

    % U3 = a1 * w2 - a2 * w1
    w2 = KernelGrid3D(r1(1),kID,   r1(2),kID_d, r1(3),kID, K);
    w1 = KernelGrid3D(r2(1),kID_d, r2(2),kID,   r2(3),kID, K);
    U(k,3) = ( sum(sum(sum(w2.*a(i1,i2,i3,1)))) - sum(sum(sum(w1.*a(j1,j2,j3,2)))) ) / h;

end