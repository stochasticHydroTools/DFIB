%-------------------------------------------------------------------------%
% spreadMAC2Dnew.m (mex version)
%
% Description: same as spreadMAC2Dnew.m, except phiweights replaced by 
% mex functions to speed up computation
% 
% By Yuanxun Bill Bao
% June 10, 2014
%-------------------------------------------------------------------------%

function f = spreadMAC2Dnew(F,X,N,h,kID,kID_d,K,L_hat)

Nx  = N(1); Ny  = N(2);
Nb  = size(X,1);
Of = zeros(Nx,Ny);
c = 1/(h*h);

for k=1:Nb
    s  = [X(k,1)/h, X(k,2)/h];
    p  = floor(s);
    r  = s-p;
    j1 = mod((p(1)-2):(p(1)+3),Nx)+1;
    j2 = mod((p(2)-2):(p(2)+3),Ny)+1;
    
    %w1 = phiweights(r(1),kID_d,K)*phiweights(r(2),kID  ,K).';
    %w2 = phiweights(r(1),kID,  K)*phiweights(r(2),kID_d,K).';
    
    % MEX function to form the 2D kernel
    w1 = KernelGrid2D(r(1), kID_d, r(2), kID  , K);
    w2 = KernelGrid2D(r(1), kID  , r(2), kID_d, K);    

    Of(j1,j2) = Of(j1,j2) + c*(F(k,2)*w1-F(k,1)*w2)/h;

    
end


%RotOf = Rot2Dnode2edge(Of,N,h);
%f(:,:,1) = PoissonSolver2D(RotOf(:,:,1),L_hat);
%f(:,:,2) = PoissonSolver2D(RotOf(:,:,2),L_hat);


LinvOf = PoissonSolver2D(Of,L_hat);
f = Rot2Dnode2edge(LinvOf,N,h);

% fix the constant
sf = (h*h)*[sum(sum(f(:,:,1))), sum(sum(f(:,:,2)))];
f0=(sum(F,1) - sf)./(h*h*Nx*Ny);

f(:,:,1) = f(:,:,1)+f0(1);
f(:,:,2) = f(:,:,2)+f0(2);



