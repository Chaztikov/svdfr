% System parameters:
Re = 2000; % Reynolds number

kzval = linspace(0.1,5,30); % spanwise wave-number
kzgrd = length(kzval);

omega = 0; % temporal frequency

y = chebfun('y');

U = 1 - y^2;    % Poiseuille flow
Uy = diff(U);
Uyy = diff(U,2);

% Looping over kz
svals = zeros(kzgrd,1);

for indz = 1:kzgrd
    A = chebop([-1 1]);
    B = chebop([-1 1]);
    C = chebop([-1 1]);

    kz = kzval(indz); k2 = kz*kz; k4 = k2*k2;

    A.op = @(x,v,eta)([1i*omega*(diff(v,2) - k2*v) - (diff(v,4)-2*k2*diff(v,2) ...
        + k4*v)/Re ;...
        1i*omega*eta + 1i*kz*Uy*v - (diff(eta,2) - k2*eta)/Re]);

    A.lbc = @(v,eta)[diff(v);v;eta];
    A.rbc = @(v,eta)[diff(v);v;eta];

    B.op = @(x,dx,dy,dz) ([ - k2*dy - 1i*kz*diff(dz);...
                            1i*kz*dx ]);
    C.op = @(x,v,eta)([ - 1i*kz*eta/k2;...
                    v ; ...
                    1i*kz*diff(v)/k2]);

    % solving for the principal singular value, use the inverse problem
    sval = svdfr(A,B,C,1,'inverse');

    % saving the largest singular value for each value of kz
    svals(indz) = 1/abs(real(sval));
    disp(indz)

end

%% Plotting the largest singular value as a function of kz at a fixed om
plot(kzval,svals,'-k');
xlabel('$k_z$');
ylabel('$\sigma_0$');
axis tight;
print('-painters','-dsvg','docs/pics/Code3_1')

%%
% System parameters:
N = 100;    % number of points for plotting
yd = chebpts(N);
clear uvec pvec
kz = 1.62; k2 = kz*kz; k4 = k2*k2;
om = 0;

A = chebop([-1 1]);
B = chebop([-1 1]);
C = chebop([-1 1]);
    
A.op = @(x,v,eta)([1i*omega*(diff(v,2) - k2*v) - (diff(v,4)-2*k2*diff(v,2) ...
        + k4*v)/Re ;...
        1i*omega*eta + 1i*kz*Uy*v - (diff(eta,2) - k2*eta)/Re]);

A.lbc = @(v,eta)[diff(v);v;eta];
A.rbc = @(v,eta)[diff(v);v;eta];

B.op = @(x,dx,dy,dz) ([ - k2*dy - 1i*kz*diff(dz);...
                        1i*kz*dx ]);
C.op = @(x,v,eta)([ - 1i*kz*eta/k2;...
                v ; ...
                1i*kz*diff(v)/k2]);

% solving for the left principal singular pair
[sfun,sval] = svdfr(A,B,C,1,'LR');

% velocities:
uvw = C(sfun(1:2,1));
u = uvw.blocks{1};
v = uvw.blocks{2};
w = uvw.blocks{3};

% Body forces:
Bad = adjointNS(B);
dxdydz = Bad(sfun(3:4,1));
dx = dxdydz.blocks{1};
dy = dxdydz.blocks{2};
dz = dxdydz.blocks{3};


% Streamfunction:
streamFun = (v)/(1i*kz);


% discretized values for plotting
pvec(:,1) = streamFun(yd); uvec(:,1) = ui(yd);
dzvec(:,1) = dz(yd);

% Getting physical fields of u and psi
zval = linspace(-2*pi/kz, 2*pi/kz, 100); % spanwise coordinate

Up = zeros(N,length(zval));   % physical value of u
Pp = zeros(N,length(zval));   % physical value of streamfunction
Dp = zeros(N,length(zval));   % physical value of body force

for indz = 1:length(zval)

    z = zval(indz);

    Up(:,indz) =  Up(:,indz) + ...
        uvec*exp(1i*kz*z) + conj(uvec)*exp(-1i*kz*z);

    Pp(:,indz) =  Pp(:,indz) + ...
        pvec*exp(1i*kz*z) + conj(pvec)*exp(-1i*kz*z);
    
    Dp(:,indz) =  Dp(:,indz) + ...
        dzvec*exp(1i*kz*z) + conj(dzvec)*exp(-1i*kz*z);

end

Up = real(Up); Pp = real(Pp); Dp = real(Dp); % only real part exist

%% Plotting the most amplified streamwise velocity structures (color plot)
pcolor(zval,yd,Up/max(max(Up))); shading interp;
colormap jet;
cb = colorbar('vert');
% Plotting the most amplified streamfunction structures (contour plot)
Ppn = Pp/max(max(Pp));
hold on
contour(zval,yd,Ppn, ...
    linspace(0.1*max(max(Ppn)),0.9*max(max(Ppn)), 4),'k--','LineWidth',1.1)
contour(zval,yd,Ppn, ...
    linspace(-0.9*max(max(Ppn)),-0.1*max(max(Ppn)), 4),'k','LineWidth',1.1)
xlabel('$z$');
ylabel('$y$');
ax = gca;
ax.XTick = -3:1:3;
print('-painters','-djpeg','docs/pics/Code3_2');
hold off;

pcolor(zval,yd,Dp/max(max(Dp))); shading interp;
colormap jet;
cb = colorbar('vert');