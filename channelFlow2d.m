% Code2: 2D channel flow example.

% Define independent variable
y = chebfun('y');

% Define the parameters of the problem:
omega = -0.3;                % frequency
Re = 2000;                   % Reynolds number
kxval = linspace(0.1,5,100); % streamwise wave-number
kxgrd = length(kxval);

U = 1 - y^2;                 % for Poiseuille flow, or U = y for Couette flow
Uy = diff(U);
Uyy = diff(U,2);

% Looping over different values of kx
svals = zeros(kxgrd,1);

parfor indx = 1:kxgrd

    kx = kxval(indx); kx2 = kx*kx; kx4 = kx2*kx2;

    % Define operator A
    A = chebop([-1 1]);
    a2 =  -( 2*kx^2/Re  +  1i*kx*U + 1i*omega );
    a0 =  kx^4/Re  +  1i*kx^3*U +  1i*kx*Uyy + 1i*omega*kx^2;
    A.op = @(y,phi) (diff(phi,4)/Re + a2*diff(phi,2) + a0*phi);

    % Specify boundary conditions
    A.lbc = @(phi) [phi;diff(phi)];
    A.rbc = @(phi) [phi;diff(phi)];

    % B operator
    B = chebop([-1 1]);
    B.op = @(y,d1,d2) (-diff(d1)  +  1i*kx*d2);

    % C operator
    C = chebop([-1 1]);
    C.op = @(y,phi) [diff(phi);-1i*kx*phi];

    % Solving for the principal singular value
    svals(indx) = svdfr(A,B,C,1);   % Compute only the largest singular value
    disp(indx);
end

%% Plot
semilogx(kxval,svals,'-k');
xlabel('$k_x$');
ylabel('$\sigma_0$');
print('-painters','-dsvg','docs/pics/Code2_1');

%%

% the OS-equation as a function of omega at kx = 1.
% System parameters:
kx = 1; % streamwise wave-number
kx2 = kx*kx; kx4 = kx2*kx2;

omval = linspace(-0.5,0,100); % temporal frequency
omgrd = length(omval);

% Looping over different values of om
svals = zeros(omgrd,1);

parfor indom = 1:omgrd
    
    omega = omval(indom);

    % Define operator A
    A = chebop([-1 1]);
    a2 =  -( 2*kx^2/Re  +  1i*kx*U + 1i*omega );
    a0 =  kx^4/Re  +  1i*kx^3*U +  1i*kx*Uyy + 1i*omega*kx^2;
    A.op = @(y,phi) (diff(phi,4)/Re + a2*diff(phi,2) + a0*phi);

    % Specify boundary conditions
    A.lbc = @(phi) [phi;diff(phi)];
    A.rbc = @(phi) [phi;diff(phi)];

    % B operator
    B = chebop([-1 1]);
    B.op = @(y,d1,d2) (-diff(d1)  +  1i*kx*d2);

    % C operator
    C = chebop([-1 1]);
    C.op = @(y,phi) [diff(phi);-1i*kx*phi];
    

    % Solving for the principal singular value
    svals(indom) = svdfr(A,B,C,1);
    disp(indom);
end

%% Plotting the largest singular value as a function of om at a fixed kx
plot(omval,svals,'-k');
xlabel('$\omega$');
ylabel('$\sigma_0$');
print('-painters','-dsvg','docs/pics/Code2_2');

%%

% System parameters:
kxval = [1 -1];        % streamwise wave-number
omval  = 0.313*[-1 1]; % temporal frequency

N = 100;               % number of Chebyshev points for plotting
yd = chebpts(N);


for n = 1:2

    omega = omval(n);
    kx = kxval(n); kx2 = kx*kx; kx4 = kx2*kx2;

    % Define operator A
    A = chebop([-1 1]);
    a2 =  -( 2*kx^2/Re  +  1i*kx*U + 1i*omega );
    a0 =  kx^4/Re  +  1i*kx^3*U +  1i*kx*Uyy + 1i*omega*kx^2;
    A.op = @(y,phi) (diff(phi,4)/Re + a2*diff(phi,2) + a0*phi);

    % Specify boundary conditions
    A.lbc = @(phi) [phi;diff(phi)];
    A.rbc = @(phi) [phi;diff(phi)];

    % B operator
    B = chebop([-1 1]);
    B.op = @(y,d1,d2) (-diff(d1)  +  1i*kx*d2);

    % C operator
    C = chebop([-1 1]);
    C.op = @(y,phi) [diff(phi);-1i*kx*phi];
    
    % Solving for the principal singular value
    [PhiAndPsi,sval] = svdfr(A,B,C,1);
    uandv = C(PhiAndPsi.blocks{1,1});   % First variable is the regular variable, phi. 
                                        % Note that C(Phi) gives the output, [u;v]
    u = uandv{1};                       % streamwise velocity
    v = uandv{2};                       % wall-normal velocity

    % discretized values for plotting
    uvec(:,n) = u(yd,1); vvec(:,n) = v(yd,1);

end

% Getting physical fields of u and v
kx = abs(kxval(1));
xval = linspace(0, 4*pi/kx, 100); % streamwise coordinate

Up = zeros(N,length(xval));   % physical value of u
Vp = zeros(N,length(xval));   % physical value of v

for indx = 1:length(xval)

    x = xval(indx);

    for n = 1:2

        kx = kxval(n);

        Up(:,indx) =  Up(:,indx) + uvec(:,n)*exp(1i*kx*x);
        Vp(:,indx) =  Vp(:,indx) + vvec(:,n)*exp(1i*kx*x);

    end

end

Up = real(Up); Vp = real(Vp); % only real part exists

%% Plotting the most amplified streamwise velocity structures
pcolor(xval,yd,Up); shading interp;
cb = colorbar('vert');
xlabel('$x$');
ylabel('$y$');
ax = gca;
ax.XTick = 0:3:12;
cb.Ticks = -3:1:3;
colormap jet;
print('-painters','-djpeg','docs/pics/Code2_3');

%% Plotting the most amplified wall-normal velocity structures
pcolor(xval,yd,Vp); shading interp;
cb = colorbar('vert');
xlabel('$x$');
ylabel('$y$');
ax = gca;
ax.XTick = 0:3:12;
cb.Ticks = -1.5:0.5:1.5;
colormap jet;
print('-painters','-djpeg','docs/pics/Code2_4');