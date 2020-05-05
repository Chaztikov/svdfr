%   Example 2: Euler Bernoulli beam example fixed at two ends
%   Code for a single value of the frequency.


% beam parameters
clear;
mu = 0.0267116; % density per unit length
EI = 128.2; % Flexural Rigidity

samp = 200;
omegas = logspace(2.5,3.5,samp);
lam = ones(10,1);
% for i = 1:samp
    omega = omegas(68);
    A = chebop([0 1]);
    B = chebop([0,1]);
    C = chebop([0,1]);

    A.op = @(x,w)(-mu*omega^2*w + EI*diff(w,4));
    B.op = @(x,q)(q);
    C.op = @(x,w)[diff(w,2);omega*w];
    A.lbc = @(w) [w;diff(w,2)];
    A.rbc = @(w) [w;diff(w,2)];

    [V,lam] = svdfr3(A,B,C,1,'LR');
    stateReg = V.blocks{1,1};
    stateReg = stateReg/sum(stateReg);
    stateAdj = -V.blocks{2,1};
    stateAdj = stateAdj/sum(stateAdj);
    outreg1 = diff(stateReg,2);
    outreg2 = omega*stateReg;
    outAdj = stateAdj;
% end
% loglog(omegas,lam);
