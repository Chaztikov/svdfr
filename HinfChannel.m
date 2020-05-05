%   Example:  computing Hinf norm
%   linearized equations for 3D channel flow.
%
%   Written by Gokul Hariharan, harih020@umn.edu

clear;
clc;
ii = sqrt(-1);
Re = 2000;
kx = 1;
kz = 1;
k2 = kx*kx + kz*kz;
k4 = k2*k2;

E = chebop([-1,1]);
A = chebop([-1,1]);
B = chebop([-1,1]);
C = chebop([-1,1]);
y = chebfun('y',[-1,1]);

U = 1 - y*y;
Uy = diff(U);
Uyy = diff(Uy);
omega = -0.4;


E.op = @(x,v,eta)([(diff(v,2) - k2*v);...
                    eta]);
A.op = @(x,v,eta)([ (diff(v,4)-2*k2*diff(v,2) ...
        + k4*v)/Re + ii*kx*Uyy*v  - ii*kx*U*(diff(v,2) - k2*v);...
        - ii*kz*Uy*v - ii*kx*U*eta + (diff(eta,2) - k2*eta)/Re]);

A.lbc = @(v,eta)[diff(v);v;eta];
A.rbc = @(v,eta)[diff(v);v;eta];

B.op = @(x,dx,dy,dz) ([-ii*kx*diff(dx) - k2*dy - ii*kz*diff(dz);...
                        ii*kz*dx - ii*kx*dz]);
C.op = @(x,v,eta)([ii*kx*diff(v)/k2 - ii*kz*eta/k2;...
                v ; ...
                ii*kz*diff(v)/k2 + ii*kx*eta/k2]);

cheboppref.setDefaults('minDimension',32);
cheboppref.setDefaults('maxDimension',1000);

% Choose one of Chebfun's discretization, ultraS works much better
cheboppref.setDefaults('discretization',@ultraS)
%cheboppref.setDefaults('discretization',@chebcolloc2)
[omega_opt,gamma_opt] = HinfNorm(E,A,B,C);

%[V,lam] = svdfr3(A,B,C,6,'LR');
