%   Example 4: frequency responses of a 2D Oldroyd-B fluid.
%   This is a code for a single value of frequency.
%   By default, this works for Poiseuille flow. For Couette flow, you may need
%   need to increase the maxDimension in cheboppref.
%   This example uses svdfri to compute the inverse problem using Matlab's
%   sparse solver, eigs.
%   Written by Gokul Hariharan, harih020@umn.edu


clear;
We = 40;
ii = sqrt(-1);
Re = 0;
beta = 0.5;
kx = 1;
k2 = kx*kx;
A = chebop([-1,1]);
B = chebop([-1,1]);
C = chebop([-1,1]);
Ctau = chebop([-1 1]);
Ctauxx = chebop([-1 1]);
Cphi = chebop([-1 1]);
y = chebfun('y',[-1,1]);

type = 'Poiseuille';
%type = 'Couette';
if strcmp(type,'Poiseuille')
        U = 1 - y*y;
    else 
        U = y;
    end
Uy = diff(U);
    Uyy = diff(Uy);
    
omega = 0;
    c = (ii * omega + 1.0 + (ii * kx * We * U));
    cy = ii * We * kx * Uy;
    cyy = ii * We * kx * Uyy;
    
    a4 = -beta + (-1 + beta)/c;
    
    a3 = (-2*(-1 + beta)*(cy - (ii)*c*kx*Uy*We))/(c^2);
    
    a2 = 2*beta*k2 - ((-1 + beta)*(-2*cy^2 - 4*kx*Uy*We*(ii*cy + kx*Uy*We) + ...
         c^2*kx*(2*kx - 3*ii*Uyy*We + 2*kx*Uy^2*We^2) + c*(cyy + 2*ii*kx*We*(Uyy + Uy*(cy + ii*kx*Uy*We)))))/...
     c^3;
    
    a1 = ((-2*ii)*(-1 + beta)*kx*(6*cy*Uy*We*(cy - ii*kx*Uy*We) + c^3*kx*Uy*We*(kx - 2*ii*Uyy*We) + ...
        c^2*(Uyy*We*(cy + 3*ii*kx*Uy*We) + kx*(ii*cy - 2*kx*(Uy^3)*We^3)) - ...
        2*c*We*(2*cy*(Uyy + ii*kx*Uy^2*We) + Uy*(cyy + kx*We*((-2*ii)*Uyy + kx*Uy^2*We)))))/c^4;
    
    a0 = (kx*(-(beta*c^4*kx^3) + 12*(-1 + beta)*cy*kx*Uy^2*We^2*(cy - ii*kx*Uy*We) + ...
        (-1 + beta)*c^3*k2*(kx - ii*Uyy*We + 2*kx*Uy^2*We^2) + ...
        (-1 + beta)*c^2*(-(cyy*kx) + ii*(-cyy + 2*k2)*Uyy*We + 2*kx*Uyy^2*We^2 + ...
           2*ii*cy*kx*Uy*We*(kx + 2*ii*Uyy*We) + 2*kx*Uy^2*We^2*(-cyy + k2 + (6*ii)*kx*Uyy*We)) + ...
        2*(-1 + beta)*c*(-2*cyy*kx*Uy^2*We^2 + cy*kx*(cy - 2*ii*kx*Uy*We)*(1 + 2*Uy^2*We^2) + ...
           ii*Uyy*We*(cy^2 + (6*ii)*cy*kx*Uy*We + 6*k2*Uy^2*We^2))))/c^4;
       a3 = a3/a4;
       a2 = a2/a4;
       a1 = a1/a4;
       a0 = a0/a4;
    b0 = 1/a4;
    b1 = -ii*kx/a4;
    A.op = @(x,v) (diff(v,4) + a3*diff(v,3) + a2*diff(v,2) + a1*diff(v) + a0*v);
    A.lbc = @(v) [diff(v);v];
    A.rbc = @(v) [diff(v);v];
    B.op = @(x,dx,dy)(b0*diff(dx) + b1*dy);
    C.op = @(x,v)[diff(v);-ii*kx*v];
    Cphi.op = @(x,v) (v);
    Ctauxx.op = @(x,v) (((4*ii*kx*We^2*Uy*Uyy)*v + (2*ii*kx+4*ii*kx*We^2*Uy^2)*diff(v) + (2*We*Uy)*diff(v,2) + (2*Uy*We*((((k2 + 2*k2*We^2*Uy^2 + ii*kx*We*Uyy)*v + diff(v,2) + Uy*We*(((2*k2*We*Uy/c)*v - (2*ii*kx/c)*diff(v))))/c))))/c);
  

cheboppref.setDefaults('maxDimension',6000);
cheboppref.setDefaults('discretization',@ultraS)
[V,lam] = svdfr3i(A,B,Ctauxx,4,'SM');
% Extract the stresses as output
stresses = Ctauxx(V.blocks{1,1});
Bad = adjointNonSquare(B);
inputs = Bad(V.blocks{2,1});
%% Plot the stressesxx
set(0,'DefaultTextInterpreter', 'latex');

plot(real(stresses),'-k',imag(stresses),'--k');
hh = legend('$\mathrm{Re}(\tau_{xx})$','$\mathrm{Im}(\tau_{xx})$','location','north');
hh.Interpreter = 'latex';
ylabel('$\tau_{xx}$','FontSize',28);
xlabel('$y$','FontSize',28);
axis tight
ax = gca;
ax.FontName = 'serif';
ax.FontSize = 28;
grid on;
ax.GridAlpha = 0.5;
grid minor;
ax.MinorGridAlpha = 0.35;
print(['stress11_We' num2str(We)],'-dpng');
ax.XLim = [0.8146    1.0000];
ax.YLim = [-105.3985   37.1812];
hh.Location = 'south';
print(['stress11_We' num2str(We) '_2'],'-dpng');
