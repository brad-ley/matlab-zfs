clearvars
close all

mu0 = 1.25663706e-6; % m kg s-2 A-2
mu_B = -9.274009994e-24; % J/T, Bohr magneton
g = 1.992; % assumed isotropic for Gd(III), as reported in Clayton et al. (2018)
hbar = 1.05457148e-34; % m2 kg / s

distance_range = linspace(1,6,201); % want 1nm to 6nm, 200 steps

sweepmid = 8.60815e3;
sweepsize = 6;
temp = 30;
bpoints = 2048;
gS = 0.00028;

data = zeros(length(distance_range), bpoints);

D_central = 1213;
D_sig = 418;
Pp_Pm = 1.6;

D_vals = linspace(-2*D_central,2*D_central,2048);
P_of_D = normpdf(D_vals,-D_central,D_sig) + ...
    Pp_Pm*normpdf(D_vals,D_central,D_sig);
PoD = P_of_D / sum(P_of_D); % normalized bimodal probability dist for D

for ii = linspace(0,length(D_vals),50)
    
Exp = struct('CenterSweep',[sweepmid sweepsize],'mwFreq',240,...
        'nPoints',bpoints,'Harmonic',0,'Temperature',temp);

% Sys = struct('g',g,'S',1/2,'gStrain',gS);
Sys = struct('g',g,'S',7/2,'D',12,'E',0);

[field, sig] = pepper(Sys, Exp);

figure
plot(linspace(-3,3,length(sig)),sig)