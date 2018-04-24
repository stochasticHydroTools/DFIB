% ----------------------------------------------------------------------- %
% Author: Yuanxun Bill Bao
% Date: July 2014
% Description: force spreading for the new divergence-free IB
% method, 3D problem
%
% Inputs:
% F - Lagrangian Forces
% X - Marker locations [x,y,z].
% N - Domain size [Nx,Ny,Nz].
% h - Cartesian grid spacing.
% kID   - kernels 'stnd3pt', 'stnd4pt', etc.
% kID_d - derivative of kernels 'stnd3pt_d', 'stnd4pt_d' , etc.
% K - Parameter for flexible kernels.
% L_hat - 3D Laplacian in Fourier space
% ----------------------------------------------------------------------- %

function f = spreadMAC3Dnew(F,X,N,h,kID,kID_d,K,L_hat)

Nx = N(1); Ny = N(2); Nz = N(3);
Nb = size(X,1);
Of = zeros(Nx,Ny,Nz,3);
LinvOf = Of;
c = 1/(h^3);        % constant in front of delta function
V = h^3*prod(N);    % vol of the box


for k = 1 : Nb
    
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
    
    
    % w2*F3 - w3*F2
    w2 = KernelGrid3D(r1(1),kID,   r1(2),kID_d, r1(3),kID,   K);
    w3 = KernelGrid3D(r1(1),kID,   r1(2),kID,   r1(3),kID_d, K);
    Of(i1,i2,i3,1) = Of(i1,i2,i3,1) + ...
                     c * (w2*F(k,3)-w3*F(k,2))/h; % divide h for derivative
    
    % w3*F1 - w1*F3
    w1 = KernelGrid3D(r2(1),kID_d, r2(2),kID,   r2(3),kID,   K);
    w3 = KernelGrid3D(r2(1),kID,   r2(2),kID,   r2(3),kID_d, K);
    Of(j1,j2,j3,2) = Of(j1,j2,j3,2) + ...
                     c * (w3*F(k,1)-w1*F(k,3))/h;
    
    % w1*F2 - w2*F1
    w1 = KernelGrid3D(r3(1),kID_d, r3(2),kID,   r3(3),kID, K);
    w2 = KernelGrid3D(r3(1),kID,   r3(2),kID_d, r3(3),kID, K);
    Of(k1,k2,k3,3) = Of(k1,k2,k3,3) + ...
                     c * (w1*F(k,2)-w2*F(k,1))/h;

end

LinvOf(:,:,:,1) = PoissonSolver3D(Of(:,:,:,1),L_hat);
LinvOf(:,:,:,2) = PoissonSolver3D(Of(:,:,:,2),L_hat);
LinvOf(:,:,:,3) = PoissonSolver3D(Of(:,:,:,3),L_hat);

f = Rot3Dedge2face(LinvOf,N,h);

% fix the constant
sf = h^3 * [ sum(sum(sum(f(:,:,:,1)))), ...
             sum(sum(sum(f(:,:,:,2)))), ...
             sum(sum(sum(f(:,:,:,3))))];
f0 = (sum(F,1) - sf) / V;

f(:,:,:,1) = f(:,:,:,1) + f0(1);
f(:,:,:,2) = f(:,:,:,2) + f0(2);
f(:,:,:,3) = f(:,:,:,3) + f0(3);


