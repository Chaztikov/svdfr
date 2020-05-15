clear;
close all;

N = 200; 
% Set parameters and base velocity
Re = 2000; % Reynolds number

y = chebfun('y');


U = 1 - y^2; % base flow (Poiseuille; U = y for Couette)
Uy = diff(U);    % derivative of U (Poiseuille; Uy = 1 for Couette)
Uyy = diff(U,2);

kxval = [1 1 -1 -1];        % streamwise wavenumber
kzval = [1 -1 1 -1];        % spanwise wavenumber
omval  = 0.385*[-1 -1 1 1]; % temporal frequency

% Looping over kx, kz, and om
uvec = zeros(N,length(kxval)); 
dzvec = zeros(N,length(kxval));

for n = 1:length(kxval)

    kx = kxval(n); kz = kzval(n); omega = omval(n);
    kx2 = kx*kx; kz2 = kz*kz;
    k2 = kx2 + kz2;

    A = chebop([-1,1]);
    B = chebop([-1,1]);
    C = chebop([-1,1]);

    A.op = @(x,v,eta)([1i*omega*(diff(v,2) - k2*v) - (diff(v,4)-2*k2*diff(v,2) ...
            + k4*v)/Re - 1i*kx*Uyy*v  + 1i*kx*U*(diff(v,2) - k2*v);...
            1i*omega*eta + 1i*kz*Uy*v + 1i*kx*U*eta - (diff(eta,2) - k2*eta)/Re]);

    A.lbc = @(v,eta)[diff(v);v;eta];
    A.rbc = @(v,eta)[diff(v);v;eta];

    B.op = @(x,dx,dy,dz) ([-1i*kx*diff(dx) - k2*dy - 1i*kz*diff(dz);...
                            1i*kz*dx - 1i*kx*dz]);
    C.op = @(x,v,eta)([1i*kx*diff(v)/k2 - 1i*kz*eta/k2;...
                    v ; ...
                    1i*kz*diff(v)/k2 + 1i*kx*eta/k2]);

    % Compute the singular function
	[sfun,sval] = svdfr3(A,B,C,1,'LR');

    % velocities
    uvw = C(sfun(1:2,:));
    ui = uvw.blocks{1}; % streamwise velocity
    vi = uvw.blocks{2}; % wall-normal velocity
    wi = uvw.blocks{3}; % spanwise velocity
    
    % Body forces:
    Bad = adjointNS(B); 
    dxdydz = Bad(sfun(3:4,:));
    dx = dxdydz.blocks{1};
    dy = dxdydz.blocks{2};
    dz = dxdydz.blocks{3};
    
    uvec(:,n) = ui(chebpts(N)); 
    dzvec(:,n) = dz(chebpts(N));
    
end

kx = abs(kxval(1)); kz = abs(kzval(1));
zval = linspace(-7.8, 7.8, 100);    % spanwise coordinate
xval = linspace(0, 12.7, 100);      % streamwise coordinate

Up = zeros(length(zval),length(xval),N);    % physical value of u
Dp = zeros(length(zval),length(xval),N);    % physical value of dz

for indx = 1:length(xval)
    x = xval(indx);
    for indz = 1:length(zval)
        z = zval(indz);
        for n = 1:4

            kx = kxval(n); kz = kzval(n);

            Up(indz,indx,:) =  squeeze(Up(indz,indx,:)) + ...
                uvec(:,n)*exp(1i*kx*x + 1i*kz*z);

            Dp(indz,indx,:) =  squeeze(Dp(indz,indx,:)) + ...
                dzvec(:,n)*exp(1i*kx*x + 1i*kz*z);
        end
    end
end



Up = real(Up); Dp = real(Dp); % only real part exist

%% Plot isosurfaces of optimal flow structures

clf
Upmax = max(max(max(Up)));
val = 0.2*Upmax;
ypt = chebpts(N);

[X,Z,Y] = meshgrid(xval,zval,ypt);
p = patch(isosurface(X,Z,Y,Up,val));
hh = isonormals(xval,zval,ypt,Up,p);
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect('auto')
view(3);
ax = gca;
xlabel('$x$');
ylabel('$z$');ax.YTick = -10:5:10;
zlabel('$y$');ax.ZTick = -1:0.5:1;
camlight

hold on
p = patch(isosurface(X,Z,Y,Up,-val));
isonormals(xval,zval,ypt,Up,p);
p.FaceColor = 'green';
p.EdgeColor = 'none';
print('-painters','-djpeg','docs/pics/Code4_2');
disp('done');
hold off
%% Plot isosurfaces of the corresponding body force
clf

dpmax = max(max(max(Dp)));
val = 0.2*dpmax;
ypt = chebpts(N);

[X,Z,Y] = meshgrid(xval,zval,ypt);
p = patch(isosurface(X,Z,Y,Dp,val));
hh = isonormals(xval,zval,ypt,Dp,p);
p.FaceColor = 'blue';
p.EdgeColor = 'none';
daspect('auto')
view(3);
ax = gca;
xlabel('$x$');
ylabel('$z$');ax.YTick = -10:5:10;
zlabel('$y$');ax.ZTick = -1:0.5:1;
camlight

hold on
p = patch(isosurface(X,Z,Y,Dp,-val));
isonormals(xval,zval,ypt,Dp,p);
p.FaceColor = 'yellow';
p.EdgeColor = 'none';
print('-painters','-djpeg','docs/pics/Code4_3');

