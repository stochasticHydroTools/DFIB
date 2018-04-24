%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% July 22, 2013
% Yuanxun Bill Bao
% 
% 2D discrete Rot operator: 
% node-centered scalar field -> edge-centered vector field
%
% Rot f = [ D2f, -D1f]
%
% For stencils: RotCurl.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = Rot2Dnode2edge(f,N,h)

Nx = N(1); Ny = N(2);

xp = [2:Nx,1]; 
yp = [2:Ny,1]; 

u = zeros(Nx,Ny,2);

u(:,:,1) =  (f(:,yp)-f)/h;

u(:,:,2) = -(f(xp,:)-f)/h;

