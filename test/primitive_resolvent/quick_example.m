%function quick_example

function [H,u,s,v] = quick_example(Re, nz, nt, beta, fund_freq, N, nsvd)

% Re = 186; % Re number
% kx = 1.5; % streamwise wavenumber
% kz = 12; % spanwise wavenumber
% cP = 10; % wave speed in plus units

% % N = 182; % Number of wall normal discretization points
% N = 10;
% nsvd = 3; % Number of singular modes from SVD

% Mitul's code allows for modification of BC's for compliant surfaces
% Leave empty for "normal" (no-slip) resolvent
yPD = [];
Atop = [];
Abot = [];

% compute mode number conversions
kz = beta*nz;
omega = fund_freq*nt;

[u,v,s,H,y,D1,D2,dy,U0] = fullChannelResolventOC(Re,0,kz,omega,N,nsvd,yPD,Atop,Abot);
% disp(LHS(1:4, 1:4))

% sigma = RIG.s(1); % first singular value

% psiu = RIG.u(1:N,1); % first singular u velocity mode
% psiv = RIG.u(1+N:2*N,1); % first singular v velocity mode
% psiw = RIG.u(1+2*N:3*N,1); % first singular w velocity mode

% phiu = RIG.v(1:N,1); % first singular u forcing mode
% phiv = RIG.v(1+N:2*N,1); % first singular v forcing mode
% phiw = RIG.v(1+2*N:3*N,1); % first singular w forcing mode

end
