function Brad_Gd_Pake(varargin)
%%% Intended to generate a Pake pattern for given input params
%%% params = (r_start, r_end, r_length, D_central, D_std, P+/P-, nMesh, T, B0, nPoints, Freq)
%%% Defaults are: (1, 6, 200, 1000, 500, 1.6, 100, 30, 8608.16, 2048, 240)
%% Begin preamble %%
format long
tic

% Debugging: Brad_Gd_Pake(1, 6, 4, 1000, 500, 1.6, 10, 30, 8608.16, 2048, 240)
PLOT = 0;
if ismac == 1
    CLUSTER = 0;
elseif ismac==0 && isunix==1
    CLUSTER = 1;
    PLOT = 0;
else
    fprintf("Not written for this OS")
end

%%% Define defaults %%%
r_start = 1; %nm
r_end = 6; %nm
r_length = 100; %200
D_central = 1000; %MHz
D_std = 500; %MHz
Pp_Pm = 1.6;
nMesh = 25; %100
T = 30; %K
B0 = 8608.16; %mT
nPoints = 2048;
Freq = 240; %GHz

%%% Define constants %%%
mu0 = 1.25663706e-6; % m kg s-2 A-2
mu_B = 9.274009994e-24; % J/T, Bohr magneton
g = 1.992; % assumed isotropic for Gd(III), as reported in Clayton et al. (2018)
hbar = 1.05457148e-34; % m2 kg / s
%% End preamble

%% Begin real function %%
varargin_reformat = zeros(length(varargin));
for ii=1:nargin %command line arguments come as strings, so we need to
    %convert to numbers in this case
    if ischar(varargin(ii))
        varargin_reformat(ii) = str2double(varargin(ii));
    elseif iscell(varargin(ii))
        varargin_reformat(ii) = varargin{ii};
    end
end

%Quick and dirty way to alter inputs from default if provided
if nargin > 0
    r_start = varargin_reformat(1);
    if nargin > 1
        r_end = varargin_reformat(2);
        if nargin > 2
            r_length = varargin_reformat(3);
            if nargin > 3
                D_central = varargin_reformat(4);
                if nargin > 4
                    D_std = varargin_reformat(5);
                    if nargin > 5
                        Pp_Pm = varargin_reformat(6);
                        if nargin > 6
                            nMesh = varargin_reformat(7);
                            if nargin > 7
                                T = varargin_reformat(8);
                                if nargin > 8
                                    B0 = varargin_reformat(9);
                                    if nargin > 9
                                        nPoints = varargin_reformat(10);
                                        if nargin > 10
                                            Freq = varargin_reformat(11);
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

distance_range = linspace(r_start,r_end,r_length);
sweeprange = 120;
coupling = -[-1/2, -1/2, 1];
data = zeros(nPoints, length(distance_range));
out = zeros(nPoints, length(distance_range));
Exp = struct('Range',[B0-sweeprange/4 B0+sweeprange/4],'mwFreq',Freq,...
        'nPoints',nPoints,'Harmonic',0,'Temperature',T); %only sweep 'up'

%%% Make probability weight functions
D_vals = linspace(-D_central-2*D_std,D_central+2*D_std,nMesh);
P_of_D = normpdf(D_vals,-D_central,D_std) + ...
    Pp_Pm*normpdf(D_vals,D_central,D_std);
PoD = P_of_D / sum(P_of_D); % normalized bimodal probability dist for D

if isempty(gcp('nocreate'))==0
    delete(gcp('nocreate'))
end

c = parcluster;
poolobj = parpool(c);
maxNumCompThreads(12);
if CLUSTER==0
    addAttachedFiles(poolobj,...
            {'/Users/Brad/Documents/MATLAB/easyspin-5.2.28/easyspin'})
elseif CLUSTER==1
    addAttachedFiles(poolobj,...
                {'/home/bdprice/matlab/easyspin'}) % for cluster run DONT
%                 FORGET TO ALSO TAKE OFF PLOTTING
end

if D_central > 0
    if PLOT == 1
        figure
        plot(D_vals, PoD)
        title('PoD')
    end
    % Outer loop
    for dd=1:length(D_vals) % loop over each D
        E_vals = linspace(0,D_vals(dd)/3,nMesh); % want same amount of separation for each value of D
        P_of_E_over_D = E_vals/D_vals(dd) - 2*E_vals.^2/D_vals(dd)^2;
        PoE = P_of_E_over_D / sum(P_of_E_over_D);
        if PLOT == 1 && dd == 1
            figure
            plot(E_vals, PoE)
            title('PoE')
        end
        for ee=1:length(E_vals) % loop over each E
            D = D_vals(dd);
            E = E_vals(ee);
            ZFS = [-1,-1,2]/3*D + [1,-1,0]*E;
            parfor ii = 1:length(distance_range)
                r = distance_range(ii) * 10^(-9);
                w_dd = mu0*mu_B^2*g^2 / (4*pi*hbar*r^3) / 1e6; % in MHz
                current_coupling = w_dd * coupling;
                Sys = struct('g',[2.00231930436153,1.992],'S',[7/2,1/2],...
                    'eeD',current_coupling,'D',[ZFS; 0,0,0]);
                [~, pake] = pepper(Sys, Exp);
                out(:, ii) = pake;
            end
            data = data + out*PoD(dd)*PoE(ee); % weight based on probability
        end
        fprintf('%.2f%% done\n', double(dd)/double(length(D_vals))*100)
    end
else
    parfor ii = 1:length(distance_range)
        r = distance_range(ii) * 10^(-9);
        w_dd = mu0*mu_B^2*g^2 / (4*pi*hbar*r^3) / 1e6; % in MHz
        current_coupling = w_dd * coupling;
        Sys = struct('g',[2.00231930436153,1.992],'S',[7/2,1/2],...
            'eeD',current_coupling);
        [~, pake] = pepper(Sys, Exp);
        data(:, ii) = pake; % normalize
    end
end

if PLOT == 1
    figure
    plot(data(:,2))
end

delete(poolobj)

savename = ['data/Pake (Raitsimring); D,stdev='...
    num2str(D_central) ',' num2str(D_std) '.txt'];

ii = 1;
newsave = savename;
while isfile(newsave)==1 % check and see if file already exists; do not overwrite if so
    erase(newsave,'.txt');
    newsave = [erase(savename,'.txt') '_' num2str(ii) '.txt'];
    ii = ii + 1;
end

if strcmp(newsave,savename)==0
    savename = newsave;
end

dlmwrite(savename,data)
    
toc
