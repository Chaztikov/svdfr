%   Example 2: Euler Bernoulli beam example fixed at two ends
%   This computes singular values for a many values of frequency, omega.
%   I've used parfor to use all cores in the computer. If you don't like
%   parallel, change parfor to for.
%   Written by Gokul Hariharan, harih020@umn.edu


% beam parameters
clear;
mu = 0.0267116; % density per unit length
EI = 128.2; % Flexural Rigidity

samp = 200;
omegas = logspace(2.5,3.5,samp);
lam = zeros(10,1);
parfor i = 1:samp
    omega = omegas(i);
    A = chebop([0 1]);
    B = chebop([0,1]);
    C = chebop([0,1]);

    A.op = @(x,w)(-mu*omega^2*w + EI*diff(w,4));
    B.op = @(x,q)(q);
    C.op = @(x,w)[diff(w,2);omega*w];
    A.lbc = @(w) [w;diff(w,2)];
    A.rbc = @(w) [w;diff(w,2)];

    lam(i) = svdfr3(A,B,C,1,'LR');
    disp(['lam = ' num2str(lam(i))]);
    disp(['omega = ' num2str(omega)]);
end
%%

loglog(omegas,lam,'-k');
ylabel('$\sigma_0$');
xlabel('$\omega$');
ax = gca;
ax.XTick = [10^2.5 10^3 10^3.5];
ax.XTickLabel = {'10^{2.5}','10^3','10^{3.5}'};
print('-painters','-dsvg','docs/pics/Euler_bernouli');
