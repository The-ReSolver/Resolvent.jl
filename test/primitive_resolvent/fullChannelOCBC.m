%% Function to impose boundary conditions and compute the resolvent
function [H0,HC] = fullChannelOCBC(LHS,RHS,N,yP,yPD,Atop,Abot)
% The NS equations have been expressed as:
% (-1i*om*M-L)[u;v;w;p] = M[fx;fy;fz;0] -> LHS [u;v;w;p] = RHS [fx;fy;fz;0]

% This function imposes the BCs corresponding to no slip and a compliant
% surface with impedance Ztop and Zbot at the top and bottom surfaces, and
% computes the resolvent

%Inputs:
% LHS: block matrix (-1i*omega*M-L)
% RHS: block matrix M
% N: grid resolution
% yP: wall-normal coordinate in plus units

% Control parameters
% yPD: detection plane distance from the wall in plus units
% Atop and Abot set the amplitude and phase of the wall-normal velocity
% ! Atop = 0 or Abot = 0 correspond to no control

%Outputs:
% H0: resolvent for no-slip
% HC: resolvent for compliant BC

% Impose no slip
LHS0 = LHS;
RHS0 = RHS;
for ni = [1,N,N+1,2*N,2*N+1,3*N]
    LHS0(ni,:)=0;
    RHS0(ni,:)=0;
    LHS0(ni,ni)=1;
end
H0 = LHS0\RHS0;

% Opposition Control BC
LHSC = LHS0;
RHSC = RHS0; 
% Find grid points closest to specified detection plane
Nbot = round(interp1(yP,1:N,min(yP)+yPD));
Ntop = round(interp1(yP,1:N,max(yP)-yPD));
% Impose blowing and suction BC: v(wall) + A*v(sensor) = 0
LHSC(N+1,N+1)=1;  LHSC(N+1,N+Ntop)=Atop;
LHSC(2*N,2*N)=1;  LHSC(2*N,N+Nbot)=Abot;
% No slip for u and w remains the same
% Compute resolvent
HC = LHSC\RHSC;
end

