% Frequency responses of the flow of 3D Oldroyd B fluid in a channel
% Written by Gokul

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
    R = 0;
    beta = 0.5;
    kx = 1;
    kz = 1;
    k2 = kx*kx + kz*kz;
    A = chebop([-1,1]);
    B = chebop([-1,1]);
    C = chebop([-1,1]);
    Ctau = chebop([-1 1]);
    Ctauxx = chebop([-1 1]);
    Cphi = chebop([-1 1]);
    y = chebfun('y',[-1,1]);

    if strcmp(type,'Poiseuille')
        U = 1 - y*y;
    else 
        U = y;
    end
    Uy = diff(U);
    Uyy = diff(Uy);
    
    T12 = Uy;
    T11 = 2*We*Uy^2;
    T12y = diff(T12);
    T11y = diff(T11);
    T11yy = diff(T11y);
    
    
    omega = 0;
    c = (ii * omega + 1.0/We + (ii * kx * U));
    cy = ii * kx * Uy;
    cyy = ii * kx * Uyy;
    
    veqD4v = -beta + (-1 + beta)/(c*We);

    veqD3v = (Complex(0,2)*(-1 + beta)*(Complex(0,1)*cy + c*kx*T12*We))/(Power(c,2)*We);

    veqD2v = (2*(-1 + beta)*(Power(cy,2) + Complex(0,2)*cy*kx*Uy + 2*Power(kx,2)*Power(Uy,2)) + Power(c,3)*(2*beta*k2 + Complex(0,1)*omega*R + Complex(0,1)*kx*R*U)*We - (-1 + beta)*Power(c,2)*(2*k2 + kx*(kx*T11 - Complex(0,3)*T12y)*We) + (-1 + beta)*c*(-cyy - Complex(0,2)*kx*Uyy + 2*kx*T12*(Complex(0,-1)*cy + kx*Uy)*We))/(Power(c,3)*We);

    veqDv = ((-1 + beta)*(-12*cy*kx*Uy*(Complex(0,1)*cy + kx*Uy) - Power(c,3)*kx*(kx*T11y + Complex(0,2)*k2*T12)*We + 2*Power(c,2)*(Power(kx,2)*((Complex(0,1)*kx*T11 + T12y)*Uy + 2*T12*Uyy)*We + cy*(Power(kx,2) + Power(kz,2) - Complex(0,1)*kx*T12y*We)) + 4*c*kx*(2*(Complex(0,1)*cy + kx*Uy)*Uyy + Complex(0,1)*Uy*(cyy + kx*T12*(Complex(0,2)*cy + kx*Uy)*We))))/(Power(c,4)*We);

    veqv = (12*(-1 + beta)*cy*Power(kx,2)*T12*Uy*(cy - Complex(0,1)*kx*Uy)*We - Power(c,4)*(k2*(beta*k2 + Complex(0,1)*omega*R) + Complex(0,1)*kx*R*(k2*U + Uyy))*We + (-1 + beta)*Power(c,3)*k2*(k2 + kx*(kx*T11 - Complex(0,1)*T12y)*We) + (-1 + beta)*Power(c,2)*(Power(kx,2)*T11y*(-cy + Complex(0,2)*kx*Uy)*We + 2*k2*kx*T12*(Complex(0,1)*cy + kx*Uy)*We - cyy*(k2 + kx*(kx*T11 + Complex(0,1)*T12y)*We) + 2*kx*Uyy*(kx*T12y*We + Complex(0,1)*(k2 + Power(kx,2)*T11*We))) + 2*(-1 + beta)*c*(Complex(0,2)*Power(kx,2)*Uy*(Complex(0,1)*cyy*T12 + kx*T12y*Uy + 2*kx*T12*Uyy)*We + Power(cy,2)*(k2 + kx*(kx*T11 + Complex(0,1)*T12y)*We) + 2*cy*kx*(-2*kx*T12*Uyy*We - Complex(0,1)*Uy*(k2 + kx*(kx*T11 - Complex(0,1)*T12y)*We))))/(Power(c,4)*We);

    etaeqD2v = (Complex(0,1)*(-1 + beta)*kz*Uy)/(Power(c,2)*We);

    etaeqDv = (Complex(0,-1)*(-1 + beta)*kz*(4*Uy*(cy - Complex(0,1)*kx*Uy) + c*(-2*Uyy + T12*(cy - Complex(0,3)*kx*Uy)*We)))/(Power(c,3)*We);

    etaeqv = (Complex(0,-1)*kz*(-(Power(c,3)*R*Uy*We) + 4*(-1 + beta)*kx*T12*Uy*(Complex(0,1)*cy + kx*Uy)*We + (-1 + beta)*c*((Power(kx,2) + Power(kz,2))*Uy + kx*T11*(Complex(0,1)*cy + 2*kx*Uy)*We - (cy*T12y + Complex(0,2)*kx*T12*Uyy)*We)))/(Power(c,3)*We);

    etaeqD2eta = -beta + (-1 + beta)/(c*We);

    etaeqDeta = (Complex(0,1)*(-1 + beta)*(Complex(0,1)*cy + kx*Uy + 2*c*kx*T12*We))/(Power(c,2)*We);

    etaeqeta = beta*k2 + Complex(0,1)*omega*R + Complex(0,1)*kx*R*U - ((-1 + beta)*kx*T12*(Complex(0,1)*cy + kx*Uy))/Power(c,2) - ((-1 + beta)*(k2 + kx*(kx*T11 - Complex(0,1)*T12y)*We))/(c*We);
    
    A.op = @(x,v,eta)([veqD4v*diff(v,4) + veqD3v*diff(v,3) + veqD2v*diff(v,2) + veqDv*diff(v) + veqv*v;...
                   etaeqD2v*diff(v,2) + etaeqDv*diff(v) + etaeqv*v + etaeqD2eta*diff(eta,2) + etaeqDeta*diff(eta) + etaeqeta*eta]);
    A.lbc = @(v,eta)[diff(v);v;eta];
    A.rbc = @(v,eta)[diff(v);v;eta];

    B.op = @(x,dx,dy,dz) ([-ii*kx*diff(dx) - k2*dy - ii*kz*diff(dz);...
                            ii*kz*dx - ii*kx*dz]);


    C.op = @(x,v,eta)([ii*kx*diff(v)/k2 - ii*kz*eta/k2;...
                    v ; ...
                    ii*kz*diff(v)/k2 + ii*kx*eta/k2]);
                
    cheboppref.setDefaults('maxDimension',1000);

    % Choose one of Chebfun's discretization, ultraS works much better
    cheboppref.setDefaults('discretization',@ultraS)
    %cheboppref.setDefaults('discretization',@chebcolloc2)
    lam = svdfr3s(A,B,C,1,'SM');
    lambdas(i) = 1/lam;
    disp(i)
    disp(1/lam)
end

%%
filename = ['samp' int2str(samp) type 'Oblique'];
save(filename);
set(0,'DefaultTextInterpreter', 'latex');

plot(Wes,abs(real(lambdas(:))),'-*k','markersize',10)
ylabel('$\sigma_0$','FontSize',28);
xlabel('$W\!e$','FontSize',28);
axis tight
ax = gca;
ax.FontName = 'serif';
ax.FontSize = 24;


plot(Wes,abs(imag(lambdas(:))),'*k','markersize',10)
ylabel('$\sigma_0$','FontSize',28);
xlabel('$W\!e$','FontSize',28);
axis tight
ax = gca;
ax.FontName = 'serif';
ax.FontSize = 24;
print('velRealObliqueP','-dpng');

function out = Complex(a,b)
    out = a + 1i*b;
end
function out = Power(a,b)
    out = a^b;
end
