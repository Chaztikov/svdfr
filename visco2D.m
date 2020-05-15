% Frequency response analysis for 2D inertialess channel flow of a viscoelastic fluid

% Set problem data
Weval = 1:20; % Weissenberg numbers
beta = 0.5;   % Viscosity ratio
omega = 0;    % temporal frequency
kx = 1;       % streamwise wavenumber
k2 = kx*kx;   % kx^2

% Flow type 
type = 'Couette'; % or 'Poiseuille';

% Singular values
svals = zeros(20,1);

for ind = 1:length(svals)
    We = Weval(ind);
    ii = sqrt(-1);
    Re = 0;
    
    A = chebop([-1,1]);  % Operator A
    B = chebop([-1,1]);  % Operator B
    C = chebop([-1,1]);  % Operator C

    y = chebfun('y');

    if strcmp(type,'Poiseuille')
        U = 1 - y*y;
    else 
        U = y;
    end
    Uy = diff(U);
    Uyy = diff(Uy);
    
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
    B.op = @(x,dx,dy)(b0*diff(dx) + b1*dy);
    C.op = @(x,v)[diff(v);-ii*kx*v];
    
    % Boundary conditions:
    A.lbc = @(v) [diff(v);v];
    A.rbc = @(v) [diff(v);v];
  
    cheboppref.setDefaults('maxDimension',2000);
    lam = svdfr(A,B,C,1,'inverse');
    svals(ind) = 1./lam';
    disp(ind)
end

%% Plot
plot(Weval,abs(real(svals(:))),'-*k')
ylabel('$\sigma_0$');
xlabel('$W\!e$');
axis tight
ax = gca;
print('-painters','-dsvg','docs/pics/Code5');
