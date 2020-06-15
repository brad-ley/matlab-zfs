clearvars

%%% Define constants %%%
mu0 = 1.25663706e-6; % m kg s-2 A-2
mu_B = 9.274009994e-24; % J/T, Bohr magneton
g = 1.992; % assumed isotropic for Gd(III), as reported in Clayton et al. (2018)
hbar = 1.05457148e-34; % m2 kg / s

g1 = 1.992;
g2 = 2.0023;
sweepmid = 8.608e3;
sweepsize = 12;
coupling = mu0*mu_B^2*g^2 / (4*pi*hbar*(2e-9)^3) / 1e6*[1/2, 1/2, -1];
temp = 30;
D = 1000;
E = 250;
ZFS = [-1,-1,2]/3*D + [1,-1,0]*E;
linewidth = [0.1 0.25];

Sys1 = struct('g',[g1, g2],'lwpp',linewidth,'S',[1/2, 7/2],'eeD',coupling,'D',[0,0,0; ZFS]);
Exp1 = struct('CenterSweep',[sweepmid sweepsize],'mwFreq',240,'nPoints',8192,...
    'Harmonic',0,'Temperature',temp);

[B1,spec1] = pepper(Sys1,Exp1);

Sys2 = struct('g',[g1, g2],'lwpp',linewidth,'S',[1/2, 7/2],'eeD',coupling,'D',[0,0,0; ZFS]);
Exp2 = struct('CenterSweep',[sweepmid sweepsize],'mwFreq',240,'nPoints',8192,...
    'Harmonic',1,'Temperature',temp);

[B2,spec2] = pepper(Sys2,Exp2);

Sys3 = struct('g',[g1, g2],'S',[1/2, 7/2],'eeD',coupling,'D',[0,0,0; ZFS]);
Exp3 = struct('CenterSweep',[sweepmid sweepsize],'mwFreq',240,'nPoints',8192,...
    'Harmonic',0,'Temperature',temp);

[B3,unbroadened] = pepper(Sys3,Exp3);

figure
plot(linspace(-sweepsize/2,sweepsize/2,length(spec1)),spec1/max(spec1))
hold on
plot(linspace(-sweepsize/2,sweepsize/2,length(spec1)),spec2/max(spec2))
legend('Spec 1', 'Spec 2')

figure
plot(linspace(-sweepsize/2,sweepsize/2,length(spec1)),unbroadened/max(unbroadened))
legend('Unbroadened')


% Sys.S = 1; Sys.g = 2; Sys.lw = 0.2;
% Sys.D = 100;
% Exp.mwFreq = 9.5; Exp.Range = [320 360]; Exp.Harmonic = 0;
% Exp.Temperature = [0.5 0.6 0.9];
% pepper(Sys,Exp);