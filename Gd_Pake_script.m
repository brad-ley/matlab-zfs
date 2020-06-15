clearvars
format long
tic

mu0 = 1.25663706e-6; % m kg s-2 A-2
mu_B = 9.274009994e-24; % J/T, Bohr magneton
g = 1.992; % assumed isotropic for Gd(III), as reported in Clayton et al. (2018)
hbar = 1.05457148e-34; % m2 kg / s

distance_range = linspace(1,6,201); % want 1nm to 6nm, 200 steps

coupling = -[-1/2, -1/2, 1];
temp = 30;
d = [0; 0];
bpoints = 2048;

data = zeros(bpoints, length(distance_range));

Exp = struct('Range',[8608.16-30 8608.16+30],'mwFreq',240,...
        'nPoints',bpoints,'Harmonic',0,'Temperature',temp);

parfor ii = 1:length(distance_range)
    
    r = distance_range(ii) * 10^(-9);
    w_dd = mu0*mu_B^2*g^2 / (4*pi*hbar*r^3) / 1e6; % in MHz
    current_coupling = w_dd * coupling;
    
    Sys = struct('g',[2.00231930436153,1.992],'S',[7/2,1/2],...
        'eeD',current_coupling);
    
    [b_field, pake] = pepper(Sys, Exp);
    
    data(:, ii) = pake/sum(pake); % normalize
end

dlmwrite('Pake Pattern (Bradneg, 1.992).txt',data)

toc

A = data;

figure(1)
plot(A(:,20))
