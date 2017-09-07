%% LOAD FILES
clc
%clear all

%This script takes velocity and pressure data from a high resolution Aquadopp
%mounted on a McLane profiler and produces profiles of dissipation rates.
%To accomplish this procedure, we employ the method outlined by Wiles et al
%(2006). 

%The subroutines (or functions) used here can be applied to any doppler
%profiler beam to obtain structure functions and dissipation rates and were
%partially inspired by Guerra & Thomson's code for 5-beam ADCPs, which is 
%rather inefficient. 
%Some of the fixes I make here should reduce computing time drastically.

%Let's load our data. Low-correlation and low-amplitude data have been
%marked already as NaNs

fpath = '/home/noel/Documents/Summer-Scripps'; % File location
%fnumber=2; %Number of files 
fname = ['testing_pitch.mat']; %File name
%We now load AQDP data as variables v1, v2, v3 and depth as variable p
load([fpath '/' fname]); 

%nbeams = 3;
<<<<<<< HEAD
blk_dist = 3; %Number of cells closest to the transducer to be removed
=======
blk_dist = 4; %Number of cells closest to the transducer to be removed
>>>>>>> f7e5b6d384b3025f8909a6a2f5aaf5d1b1418ed0
aqdp = rmfield(aqdp,{'v3','v1'}); 
%Let's load velocity data into the GPU to speed up computing
%aqdp.v1 = gpuArray(aqdp.v1(:,blk_dist+1:end)); 
%vstar = gpuArray(aqdp.v2(:,blk_dist+1:end)); 
%aqdp.v3 = gpuArray(aqdp.v3(:,blk_dist+1:end)); 
vstar = aqdp.v2(:,blk_dist+1:end);
aqdp = rmfield(aqdp, {'v2'});

<<<<<<< HEAD
%% REDUCE MEMORY LOAD

distances = (0:33)*0.022 + 0.167; %Distance between bins and transducer
distances = distances(blk_dist+1:end); 

%Below we define the number of profiles over which the products between 
%velocity differences are going to be averaged when computing structure
%functions. These sets of profiles are called ensembles. 
ens_length = 15; 
=======
%------------------------------------------------ COMMENTING SECTION

distances = (0:33)*0.022 + 0.167; %Distance between bins and transducer
distances = distances(blk_dist+1:end); 

%Below we define the number of profiles over which the products between 
%velocity differences are going to be averaged when computing structure
%functions. These sets of profiles are called ensembles. 
ens_length = 5; 
>>>>>>> f7e5b6d384b3025f8909a6a2f5aaf5d1b1418ed0

%The following array will tell us where all ensembles begin and will also
%begin to define indexing for structure function and dissipation data

starts = 1:ens_length:length(aqdp.p); 
starts = starts(1:end-150);

%Before computing dissipation at each one of these moments, let's see what
%their place (depth, time) will be:
diss_time = aqdp.yday(starts); 
diss_depth = aqdp.p(starts); %Lower bound of profile

%Set start and end date as day number for computation:
<<<<<<< HEAD
start_dt = 32.75;
end_dt = 32.7581;
%end_dt = 34.75;

=======
start_dt = 34.035;
end_dt = 35.55;
>>>>>>> f7e5b6d384b3025f8909a6a2f5aaf5d1b1418ed0
[c start_dt] = min(abs(diss_time-start_dt));
[c end_dt] = min(abs(diss_time - end_dt)); 

diss_time = diss_time(start_dt:end_dt); 
diss_depth = diss_depth(start_dt:end_dt); 
starts = starts(start_dt: end_dt); 

%Let's obtain structure functions at many times and store them all
%together at Deez to make fits and derive dissipation later:
%n_profiles = 1e5; %number of dissipation estimates
n_profiles = length(starts)
<<<<<<< HEAD

r_min = 1; r_max = 30;

%Making sure everything will work correctly:
if ge(r_min, r_max)
    warning('r_min cannot be equal or greater than r_max')
elseif r_max > size(vstar, 2)- r_min
    warning('r_max is greater than what may be computed by your selection of r_min and the size of your velocity matrix');
end

%% STRUCTURE FUNCTION ESTIMATION

Deez = NaN(n_profiles,r_max - r_min + 1); 

tic

for time = 1:n_profiles
    Deez(time,:) = AQHR_StructureFunction(vstar(starts(time): ...
        starts(time)+ens_length-1,:), r_min, r_max, size(vstar,2), 0.4,); 
=======

r_min = 2; r_max = 23;

%Making sure everything will work correctly:
if ge(r_min, r_max)
    warning('r_min cannot be equal or greater than r_max')
elseif r_max > size(vstar, 2)- r_min
    warning('r_max is greater than what may be computed by your selection of r_min and the size of your velocity matrix');
>>>>>>> f7e5b6d384b3025f8909a6a2f5aaf5d1b1418ed0
end

<<<<<<< HEAD
k = 'done'
dissipation = NaN(n_profiles,1);
noise = NaN(n_profiles,1); 
tic
radius = ((r_min:r_max)*0.022);

for structures = 1:n_profiles
    %[fit, stats] = robustfit(radius, Deez(structures,:));
    fit = robustfit(radius.^(2/3), Deez(structures,:),'ols'); 
    %fit = regress(Deez(structures,:), [ones(size(radius)), radius.^(2/3)]);
    dissipation(structures) = (abs(fit(2))/2.1)^(3/2);
    noise(structures) = fit(1); 
    %sigma(structures) = stats.s;
    %structures
=======
%---------------------------------------------- STRUCTURE FUNCTION

% Deez = NaN(n_profiles,r_max - r_min + 1); 
% 
% tic
% %parpool(2)
% %Matrices = NaN(12, r_max - r_min + 1,n_profiles);
% 
% for time = 1:n_profiles
%     [Deez(time,:)] = AQHR_StructureFunction(vstar(starts(time): ...
%         starts(time)+ens_length-1,:), r_min, r_max); 
%     %Deez(time,:) = AQHR_StructureFunction(aqdp.v2(starts(time): ...
%     %    starts(time)+ens_length-1,:))
%     %time
% end
% toc
% 
% k = 'done'
% dissipation = NaN(n_profiles,1);
% noise = NaN(n_profiles,1); 
% %sigma = NaN(n_profiles,1);
% %n_profiles = length(Deez);
% tic
% radius = ((r_min:r_max)*0.022);
% 
% for structures = 1:n_profiles
%     %[fit, stats] = robustfit(radius, Deez(structures,:));
%     fit = robustfit(radius.^(2/3), Deez(structures,:),'ols'); 
%     %fit = regress(Deez(structures,:), [ones(size(radius)), radius.^(2/3)]);
%     dissipation(structures) = (abs(fit(2))/2.1)^(3/2);
%     noise(structures) = fit(1); 
%     %sigma(structures) = stats.s;
%     %structures
% end
% 
% toc

%------------------------------------------------- SPECTRAL METHOD

bindistance = 0.022;
Corrz = NaN(n_profiles, r_max-r_min +1);
meanvs = NaN(n_profiles,1); 

%Compute correlation functions for all profiles at times starts and of
%length ens_length:
tic
for time = 1:n_profiles
    [Corrz(time,:), meanvs(time)] = AQHR_correlation(vstar(starts(time): ... 
        starts(time)+ ens_length-1,:), r_min, r_max, 1);
>>>>>>> f7e5b6d384b3025f8909a6a2f5aaf5d1b1418ed0
end
toc
%%
%------------------------------------------------- SPECTRAL METHOD

<<<<<<< HEAD
bindistance = 0.022;

window_length = 19; %number of bins taken for spectral estimation

%Spectra = NaN(n_profiles, floor((r_max - r_min +1)/2-1)); %Array to save spectra
Spectra = NaN(n_profiles, floor(window_length/2));

%Compute correlation functions for all profiles at all times in an ensemble
%and average the spectra for all those moments. 
tic
for time = 1:n_profiles
    %time
    Spectra(time,:) = AQHR_correlation(vstar(starts(time):starts(time)+ ...
        ens_length-1,:),r_min,r_max, 1,window_length, 0.4,bindistance);
    
    %AQHR_correlation computes the average PSD of correlation functions in
    %an ensemble given by the range of vstar.
end
toc

dissipation = NaN(n_profiles, 1); 
S = NaN(n_profiles, 2);
%wavenums = (((2*pi)./(bindistance*(r_min:r_min+floor(window_length/2)-1)))).^(-5/3);
%wavenums = (1./(bindistance*(r_min:r_max))).^(-5/3);
wavenums = ((1:size(Spectra,2))/(window_length*bindistance));
wavefit = wavenums.^(-5/3);

tic
%We have to do the fit only for a smaller portion of the data, so we pick
%where to make the cut: 
cut_low = 6; cut_hi = 14;

for correlation = 1:n_profiles
    S(correlation,:) = robustfit(wavefit, Spectra(correlation,:), 'ols');  
   
    %Negative slopes will give us imaginary dissipations, so we'll get rid
    %of those next
end
toc
=======

radius = ((r_min:r_max)*bindistance);
dissipation = NaN(n_profiles, 1); 
espectros = NaN(n_profiles, 9); 
%Now compute the power spectra of all those correlation functions
tic

wavenums = fliplr(log10((2*pi)./(bindistance*(r_min:r_min+8)))).^(-5/3);

for correlation = 1:n_profiles
    spectrum = obmPSpec(Corrz(correlation,:), bindistance, 19, 0.5);
    spectrum = spectrum.psd*(meanvs(correlation)/sum(spectrum.psd)); 
    
    espectros(correlation,:) = spectrum;
    fit = robustfit(wavenums, ...
        spectrum); 
    dissipation(correlation) = abs(fit(2));
end
toc

dissipation2 = (dissipation.*(55/18)).^(3/2)*(9*0.4/8);
>>>>>>> f7e5b6d384b3025f8909a6a2f5aaf5d1b1418ed0

%dissipation(dissipation < 0) = NaN; 
dissipation = (9*0.4/8).*((55/18)*abs(S(:,2))).^1.5;

