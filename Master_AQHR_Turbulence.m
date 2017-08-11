clc
clear all

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
fname = ['raw_alongbeam.mat']; %File name
%We now load AQDP data as variables v1, v2, v3 and depth as variable p
load([fpath '/' fname]); 

%nbeams = 3;
%deltaZ = 0.03; %cell height in meters
aqdp = rmfield(aqdp,{'v1','v3'}); 

%Let's load velocity data into the GPU to speed up computing
%We blank the first 3 bins so load starting on the 4th.
%aqdp.v1 = gpuArray(aqdp.v1(:,5:end)); 
aqdp.v2 = gpuArray(aqdp.v2(:,5:end)); 
%aqdp.v3 = gpuArray(aqdp.v3(:,5:end)); 
%aqdp.p = gpuArray(aqdp.p);

%Below we define the length of ensembles over which the products between 
%velocity differences are going to be averaged when computing structure
%functions. 
ens_length = 12; %Each ensemble will include 6 independent profiles over 
%which height also varies, so we need to change it

%The following array will tell us where all ensembles begin and will also
%begin to define indexing for structure function and dissipation data

starts = 1:ens_length:length(aqdp.p); 
starts = starts(1:end-100);

%Before computing dissipation at each one of these moments, let's see what
%their place (depth, time) will be:
diss_time = aqdp.yday(starts); 
diss_depth = aqdp.p(starts); %Lower bound of profile

%Let's obtain structure functions at many times and store them all
%together at Deez to make fits and derive dissipation later:
%n_profiles = 1e5; %number of dissipation estimates
n_profiles = length(starts);
Deez = NaN(n_profiles,27); 

%starting_point = floor((length(starts))/2 - n_profiles/2)-100;
tic

for time = 1:n_profiles
    Deez(time,:) = gather(AQHR_StructureFunction(aqdp.v2(starts(time): ...
        starts(time)+ens_length-1,:), 0.02)); 
    %time
end
toc

k = 'done'
dissipation = NaN(n_profiles,1);
noise = NaN(n_profiles,1); 
sigma = NaN(n_profiles,1);
n_profiles = length(Deez);
tic
radius = ((1:27)*0.02+0.04).^(2/3);

for structures = 1:n_profiles
    %[fit, stats] = robustfit(radius, Deez(structures,:));
    fit = robustfit(radius, Deez(structures,:)); 
    dissipation(structures) = (abs(fit(2))/2.1)^(3/2);
    noise(structures) = abs(fit(1)); 
    %sigma(structures) = stats.s;
    %structures
end

toc

%diss_depth = diss_depth(1:end-100); 
%diss_time = diss_time(1:end); 

save('newresult6.mat', 'Deez','dissipation','noise', ...
    'diss_depth','diss_time'); 




