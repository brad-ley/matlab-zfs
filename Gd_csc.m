function gd_csc(varargin)
%%% Intended to generate a Pake pattern for given input params
%%% params = (r_start, r_end, r_length, D_central, D_std, P+/P-, nMesh, T, B0, nPoints, Freq)
%%% Defaults are: (1, 6, 200, 1000, 500, 1.6, 100, 30, 8608.16, 2048, 240)
%% Begin preamble %%
format long
fprintf('Begin at: %s\n', datestr(now))
tStart = tic;

% Debugging: gd_csc(1, 6, 4, 1000, 500, 1.6, 10, 30, 8608.16, 2048, 240)
% Real test: gd_csc(1,6,200,d,dstd,pppm,100,30,8608.16,2048,240)
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
r_length = 200;
D_central = 1000; %MHz
D_std = 500; %MHz
Pp_Pm = 1.6;
nMesh = 100;
T = 30; %K
% B0 = 8608.16; %mT -- better to find B0 using constants and input frequency
nPoints = 2048;
Freq = 240; %GHz

%%% Define constants %%%
mu0 = 1.25663706e-6; % m kg s-2 A-2
mu_B = 9.274009994e-24; % J/T, Bohr magneton
g = 1.992; % assumed isotropic for Gd(III), as reported in Clayton et al. (2018)
hbar = 1.05457148e-34; % m2 kg / s

B0 = hbar*2*pi*Freq*1e12/(g*mu_B);
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

sweeprange = 120;
coupling = [1/2, 1/2, -1]; %eeDz negative because spins are most likely to align antiparallel

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
                {'/home/bdprice/matlab/easyspin'})
end

if abs(D_central) ~= 0
    %% Make probability weight functions
    D_vals = linspace(-D_central-2*D_std,D_central+2*D_std,nMesh);
    P_of_D = normpdf(D_vals,-D_central,D_std) + ...
        Pp_Pm*normpdf(D_vals,D_central,D_std);
    PoD = P_of_D / sum(P_of_D); % normalized bimodal probability dist for D
    
    if PLOT == 1
        figure
        plot(D_vals, PoD)
        title('PoD')
    end
    %% Outer loop
    d_out = zeros(nMesh, nPoints, r_length);
    e_spacing = max(D_vals)/3/(nMesh/2);
    parfor dd=1:nMesh % loop over each D
       out = zeros(nPoints, r_length);
       Exp = struct('Range',[B0-sweeprange/4 B0+sweeprange/4],'mwFreq',Freq,...
       'nPoints',nPoints,'Harmonic',0,'Temperature',T); %only sweep 'up'
       if D_vals(dd) ~= 0
           E_vals = 0:sign(D_vals(dd))*e_spacing:D_vals(dd)/3;
           runrest = 1;
       else
           E_vals = 0;
           runrest = 0;
       end
       if runrest~=0
           for ee=1:length(E_vals) % loop over each E
%          if D_vals(dd) == 0
%              E_vals = 0; % only defined so MatLab doesn't warn me about clearing it in every iteration of parfor loop
               P_of_E_over_D = E_vals/D_vals(dd) - 2*E_vals.^2/D_vals(dd)^2;
               if sum(P_of_E_over_D) ~= 0
                   PoE = P_of_E_over_D / sum(P_of_E_over_D);
               else
                   PoE = 0;
               end
               if PLOT == 1 && dd == 1
                   figure
                   plot(E_vals, PoE)
                   title('PoE')
               end
               
               distance_range = linspace(r_start,r_end,r_length);

               for ii=1:r_length
                   D = D_vals(dd);
                   E = E_vals(ee);
                   ZFS = [-1,-1,2]/3*D + [1,-1,0]*E;
                   r = distance_range(ii) * 10^(-9);
                   w_dd = mu0*mu_B^2*g^2 / (4*pi*hbar*r^3) / 1e6; % in MHz
                   current_coupling = w_dd * coupling;
                   Sys = struct('g',[1.992,1.992],'S',[7/2,7/2],...
                       'eeD',current_coupling,'D',[ZFS; 0,0,0]);
                   [~, pake] = pepper(Sys, Exp);
                   out(:,ii) = out(:,ii) + reshape(pake*PoD(dd)*PoE(ee),2048,[]);
               end
               fprintf('%.1f%% done d=%d of %d\n', [100*ee/length(E_vals), dd, nMesh])
           end
       end
       d_out(dd, :, :) = out;
       fprintf('approx %.1f%% done\n', 100*dd/nMesh)
    end
    data = reshape(sum(d_out, 1),[nPoints, r_length]);
else
    data = zeros(nPoints, r_length);
    Exp = struct('Range',[B0-sweeprange/4 B0+sweeprange/4],'mwFreq',Freq,...
        'nPoints',nPoints,'Harmonic',0,'Temperature',T); %only sweep 'up'
    parfor ii = 1:r_length
        distance_range = linspace(r_start,r_end,r_length);
        r = distance_range(ii) * 10^(-9);
        w_dd = mu0*mu_B^2*g^2 / (4*pi*hbar*r^3) / 1e6; % in MHz
        current_coupling = w_dd * coupling;
        Sys = struct('g',[2.00231930436153,1.992],'S',[7/2,1/2],...
            'eeD',current_coupling);
        [~, pake] = pepper(Sys, Exp);
        data(:, ii) = pake; % normalize
    end
end

B = transpose(linspace(B0-sweeprange/4, B0+sweeprange/4, nPoints));

for ii=1:length(data(1,:))
    data(:,ii)=data(:,ii)/trapz(B,data(:,ii)); %normalize area to 1 over the field sweep axis
end

if PLOT == 1
    figure
    plot(data(:,2))
end

delete(poolobj)

savename = ['data/D,stdev='...
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

fprintf('Begin at: %s\n', datestr(now))
tEnd = toc(tStart);
fprintf('%d minutes and %.1f seconds\n', floor(tEnd/60), rem(tEnd,60));
dlmwrite(savename,data)
