%% Primitive-variable resolvent code for channel flow
% Written by Mitul Luhar 01/06/2015
% For details on the primitive-variable resolvent analysis see:
% Luhar, Sharma, McKeon (2014), Opposition Control within the Resolvent 
% Analysis Framework, J. Fluid Mech., 749:597-626
% This particular code performs phase-shifted opposition control on both
% walls

% NOTE: Eventually, modularize each of the cells below

% ! ask Sean both these questions tomorrow:
% !     - does it matter if the last line of L (continuity) is positive or negative???

% function [RIG,CON,H0,y,D1,D2,dy,U0] = fullChannelResolventOC(Re,kx,kz,omega,N,nsvd,yPD,Atop,Abot)
function [u0,v0,s0W,H0,y,D1,D2,dy,U0] = fullChannelResolventOC(Re,kx,kz,omega,N,nsvd,yPD,Atop,Abot)
%% Inputs
% N:    Grid resolution
% nsvd: Number of singular modes to compute
% 
% % Flow Parameters
% Re:   Reynolds number (Re = u_tau h / nu)
% kx:   streamwise wavenumber (normalized by half-height,h)
% kz:   spanwise wavenumber (normalized by half-height,h)
% cP:   Wave speed in plus units (normalized by u_tau)
% omega = cP*kx;  % Radian frequency
% cP = omega/kx;

% Control parameters
% yPD: detection plane distance from the wall in plus units
% Atop and Abot set the amplitude and phase of the wall-normal velocity
% Atop = 0 or Abot = 0 correspond to no control

%% Set up coordinate system and estimate mean velocity profile
[y, DM] = chebdif(N, 2);
[~, dy] = clencurt(N - 1);
D1 = DM(:, :, 1);
D2 = DM(:, :, 2);
U0 = y;
% [y,D1,D2,dy,U0] = channelMeanVel(Re,N);
% y:[-1 1]
% D1: First differential
% D2: Second differential
% dy: Integration weights for scaling the resolvent
% U0: Mean velocity profile

%% Ok, now set up the linear NS operator and eventually the resolvent
% A few basic matrices
I   = eye(N);
Z   = zeros(N);

% Calculate mean shear
dU0= D1*U0;
U0 = diag(U0);
dU0 = diag(dU0);

% Create important block components
ikU0 = 1i*kx*U0;
LAP  = -(kx^2)*I - (kz^2)*I + D2;

% Block matrix L representing linearized NS equations
% The last column is the pressure gradient, last row is continuity
L1 = [-ikU0+LAP/Re, -dU0        , Z          , -1i*kx*I]; %u (ax.)
L2 = [Z           , -ikU0+LAP/Re, Z          ,      -D1]; %v (rad.)
L3 = [Z           , Z           ,-ikU0+LAP/Re, -1i*kz*I]; %w (az.)
L4 = [-1i*kx*I    ,-D1          ,-1i*kz*I    ,        Z]; %contin.
L  = [L1; L2; L3; L4];

% Mass Matrix
M = [I Z Z Z; Z I Z Z; Z Z I Z; Z Z Z Z];

% The governing equation reads: (-i*om*M - L) [u;v;w;p] = M [fx;fz;fw;0]
% ! the original code had the time derivative negated, I need to figure out why!!!
LHS = 1i*omega*M-L;
RHS = M;

%% Apply boundary conditions and compute resolvent
yP = Re*y;
[H0,HC] = fullChannelOCBC(LHS,RHS,N,yP,yPD,Atop,Abot);

% H0: basic resolvent
% HC: resolvent accounting for opposition control BC

% Scale resolvent and perform SVD
IW   = sqrtm(diag(dy));
iIW  = I/IW;
sqW  = [ IW Z Z Z; Z  IW Z Z; Z Z  IW Z; Z Z Z Z];
isqW = [iIW Z Z Z; Z iIW Z Z; Z Z iIW Z; Z Z Z Z];

% Weighted resolvents
H0W = sqW*H0*isqW;
HCW = sqW*HC*isqW;

% Singular value decomposition
[u0W,s0W,v0W] = svds(H0W,nsvd); 
u0 = isqW*u0W;
v0 = isqW*v0W;
s0 = diag(s0W);

% Because of the l2 norm used to scale the resolvent, we do not have any
% pressure data.  Calculate pressure modes using the un-scaled resolvent.
u0 = H0*v0;
u0 = u0*diag(1./s0);
RIG.u = u0;
RIG.s = s0;
RIG.v = v0;
end
