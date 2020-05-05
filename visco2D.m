%   Example 4: frequency responses of a 2D Oldroyd-B fluid.
%   By default, this works for Poiseuille flow. For Couette flow, you may need
%   need to increase the maxDimension in cheboppref.
%   This example uses svdfri to compute the inverse problem using Matlab's
%   sparse solver, eigs.
%   Written by Gokul Hariharan, harih020@umn.edu



%% Set up data
samp = 20;
Wes = zeros(samp,1);
for i = 1:samp
    Wes(i) = i;
end
%type = 'Poiseuille';
type = 'Couette';
lambdas = zeros(samp,1);
lambdasT = zeros(samp,1);

parfor i = 1:samp
    We = Wes(i);
    ii = sqrt(-1);
    Re = 0;
    beta = 0.5;
    kx = 1;
    k2 = kx*kx
    A = chebop([-1,1]);
    B = chebop([-1,1]);
    C = chebop([-1,1]);

    y = chebfun('y',[-1,1]);

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
  
    cheboppref.setDefaults('maxDimension',2000);
    cheboppref.setDefaults('discretization',@ultraS)
    lam = svdfr3i(A,B,C,1,'SM');
    lambdas(i) = 1./lam';
    disp(i)
end

%% Plot
plot(Wes,abs(real(lambdas(:))),'-*k')
ylabel('$\sigma_0$');
xlabel('$W\!e$');
axis tight
ax = gca;
print('-painters','-dsvg','docs/pics/velP');
