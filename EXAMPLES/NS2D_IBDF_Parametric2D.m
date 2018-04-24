clear; close all; restoredefaultpath;
addpath('../Kernels_MEX',...
        '../DFIB_SpreadInterp2D_MEX',...
        '../DFIBsolver2D');

L=5; 
mu=0.15;
rho = 1;

% membrane parameters
alpha = 1/5;
R = alpha * L;
epsilon = 0.05;
p = 2;

MembranePara.omega0=10;
MembranePara.Kc=10;
MembranePara.tau=0.5;
MembranePara.R = R;


Kernel = {'flex6pt','flex6pt_d', (59/60)*(1-sqrt(261/3481))};
% Kernel = {'stnd4pt','stnd4pt_d', []};
% Kernel = {'stnd3pt','stnd3pt_d', []};
% Kernel = {'bspline4pt', 'bspline4pt_d',[]};
% Kernel = {'new4ptC2','new4ptC2_d', []};
% Kernel = {'bspline6pt','bspline6pt_d', []};
% Kernel = {'flex5pt', 'flex5pt_d', (38 - sqrt(69))/60};


showplot = 'on';

% Eulerian grid
Nx=128; Ny=Nx;
N=[Nx,Ny];
h=L/Nx;

% % Lagrangian grid 
% Nb=6*Nx;
% ds=2*pi/Nb;
% s =(0:Nb-1)*ds;
% %alpha = 5/28; beta = 7/20;
% %X = [alpha*cos(s'), beta*sin(s')]*L+L/2;
% alpha = 1/5; beta = 1/5;
% epsilon=0.05;
% X = [alpha*cos(s').*(1+epsilon*cos(p*s')), beta*(1+epsilon*cos(p*s')).*sin(s')]*L+L/2;

% Lagrangian grid 
% Nb=3*Nx;
% ds=2*pi/Nb;
% s =(0:Nb-1)*ds;
%alpha = 5/28; beta = 7/20;
%X = [alpha*cos(s'), beta*sin(s')]*L+L/2;

MpC = 2; % number markers per cell
Ns = round( MpC*(2*pi*R)/(L/Nx) );
ds=2*pi/Ns;
s =(0:Ns-1)*ds;
X = [alpha*cos(s').*(1+epsilon*cos(p*s')), alpha*(1+epsilon*cos(p*s')).*sin(s')]*L+L/2; 

% initial velocity
u=zeros(Nx,Ny,2);

% time step
tend = 25;
dt   = h/20;
Nt   = floor(tend/dt);
dt   = tend/Nt;
tt   = 0:dt:tend;
Nf   = Nt/8;


tic
%[uIBDF,XIBDF,areaIBDF,amplIBDF,areaTracer]=NS2D_IBDFparametric_midpoint(L,N,mu,Kernel,X,u,tend,dt,Nf,showplot);
[uIBDF,XIBDF,areaIBDF]=Timestepping_Parametric(L,N,mu,rho,MembranePara,Kernel,X,u,tend,dt,Nf,showplot);
toc
