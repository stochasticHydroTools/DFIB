clear; close all; restoredefaultpath;
addpath('../Kernels_MEX',...
        '../DFIB_SpreadInterp2D_MEX',...
        '../DFIBsolver2D');
    
profile off;    

L=1; % domain: [0,L]x[0,L]
mu=0.1; % fluid viscosity
rho=1; % density

% Choice of IB kernel
Kernel = {'flex6pt','flex6pt_d', (59/60)*(1-sqrt(261/3481))};
%Kernel = {'stnd4pt','stnd4pt_d', []};
%Kernel={'bspline4pt', 'bspline4pt_d',[]};
%Kernel={'bspline6pt', 'bspline6pt_d',[]};
%Kernel = {'flex5pt', 'flex5pt_d', (38 - sqrt(69))/60};

showplot = 'on';

% Eulerian grid
Nx=128; Ny=Nx;
N=[Nx,Ny];
h=L/Nx;

% Lagrangian grid
alpha=1/4; beta=1/4;
MpC = 4; % number markers per cell
Ns = round(MpC*(2*pi*alpha)/L*Nx); % number of Lagrangian points
ds=2*pi/Ns;
s =(0:Ns-1)*ds;
X0 = [alpha*cos(s'), beta*sin(s')]*L+L/2;

% initial velocity
u=zeros(Nx,Ny,2);

% time step
tend = 1;
dt   = h;
Nt   = floor(tend/dt);
dt   = tend/Nt;
tt   = 0:dt:tend;
Nf   = Nt/8; % number of frames to plot

tic;


[uIBDF,XIBDF,areaIBDF] = Timestepping_Circle(L,N,mu,rho,Kernel,X0,u,tend,dt,Nf,MpC,showplot);
%[areaIBDF,areaTracer,uT,XT,XtracerT,TH,Fpol]=NS2D_IBDF_midpoint(L,N,mu,Kernel,X0,u,tend,dt,Nf,MpC,showplot);
% save(['NS2D-IBDF-circle-spacing-',num2str(MpC),'-',Kernel{1},'-midpoint.mat'], ...
%       'L','mu','Nx','MpC','Ns','ds','s',...
%       'alpha','beta','X0','tend','dt','tt', 'Kernel')
% save(['NS2D-IBDF-circle-spacing-',num2str(MpC),'-',Kernel{1},'-midpoint.mat'], ...
%       'areaIBDF','areaTracer','uT','XT','XtracerT','TH','Fpol','-append');


% [areaIBDF,areaTracer,areaTracer_para,TH,F,Xtracer]=NS2D_IBDF_AB2(L,N,mu,Kernel,X,u,tend,dt,Nf,MpC,showplot);
% save(['NS2D-IBDF-circle-spacing-',num2str(MpC),'-',Kernel{1},'-AB2-dt2-morepts0.mat'],'TH','F','areaIBDF','areaTracer','tt','Xtracer','areaTracer_para');
% 

toc

% profile viewer
% p = profile('info');
% profsave(p,'profile_results')
% profile off;
