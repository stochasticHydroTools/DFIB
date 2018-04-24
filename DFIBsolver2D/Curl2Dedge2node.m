%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% July 22, 2013
% Yuanxun Bill Bao
% 
% 2D discrete Curl operator in 2D: 
% edge centered vector field -> node-centered scalar field
%
% Curl u = D1u2 - D2u1
%
% For stencils: RotCurl.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = Curl2Dedge2node(u,N,h)

Nx = N(1); Ny = N(2);

xm = [Nx,1:Nx-1];
ym = [Ny,1:Ny-1];

f = (u(:,:,2)-u(xm,:,2))/h - (u(:,:,1)-u(:,ym,1))/h;
