%% LOAD FILES
%Ultimate script to estimate dissipation from AQDP profiles through
%structure function and spectral techniques.

clc
clear all

fpath = '/home/noel/Documents/Summer-Scripps'; % File location 
fname = ['testing_pitch.mat']; %File name

%We now load AQDP data as variables v1, v2, v3 and depth as variable p
load([fpath '/' fname]); 


blk_dist = 2; %Number of cells closest to the transducer to be removed
aqdp = rmfield(aqdp,{'v3','v1'}); 
vstar = aqdp.v2(:,blk_dist+1:end);
aqdp = rmfield(aqdp, {'v2'});

%% SELECT ENSEMBLES AND TIMES TO AVERAGE

%Set the number of profiles included to compute only one spectrum:
ens_length = 15;

%Indices of data spaced by ens_length: 
starts = 1:ens_length:length(aqdp.p); 
starts = starts(1:end-1);

diss_time = aqdp.yday(starts);
diss_depth = aqdp.p(starts);

%Set start and end date as day number for computation:
start_dt = 32.75;
end_dt = 32.7581;

[c start_dt] = min(abs(diss_time-start_dt));
[c end_dt] = min(abs(diss_time - end_dt)); 

diss_time = diss_time(start_dt:end_dt); 
diss_depth = diss_depth(start_dt:end_dt); 
starts = starts(start_dt: end_dt); 
n_profiles = length(starts);
%We are now left with time and pressure data corresponding to each of the
%ensembles for which we're going to estimate epsilon.

%% SPECTRAL METHOD
% Spectra are obtained for all ensembles and for all wavenumbers possible. 
% An ensemble is a set of consecutive AQDP profiles.

spectra = NaN(n_profiles, floor(size(vstar,2)/2-1));
r_min = 1;
r_max = size(vstar,2);

wavenums = (1:floor(size(vstar,2)/2-1))/(0.022*size(vstar,2));
wavefit = wavenums.^(5/3); wavefit2 = wavefit.^(-1);

for ensemble = 1:n_profiles
    spectra(ensemble,:) = norm_spectrum(vstar(starts(ensemble): ...
        starts(ensemble)+ens_length-1,:),1);  
    spectra(ensemble,:) = spectra(ensemble,:)./0.52;
end
 
S = NaN(n_profiles,2); 
Saux = S;

for correlation = 1:n_profiles  
    S(correlation,:) = robustfit(wavefit2, spectra(correlation,:), 'ols',1,'off'); 
    Saux(correlation,:) = robustfit(wavefit2, spectra(correlation,:)*0.52,'ols',1,'off');
end
k = 5
dissipation2 = (9*0.4/8).*((55/18)*abs(Saux(:,2))).^1.5;
dissipation = S(:,2).^(1.5); 

k = sgolayfilt(log10(dissipation), 2, 35); 
k2 = sgolayfilt(log10(dissipation2), 2, 35); 


plot(diss_time, k,'b');
hold on
plot(diss_time, k2, 'r'); 
grid on

%% STRUCTURE FUNCTION 

r_min = 2; r_max = 26;

Deez = NaN(n_profiles, r_max - r_min + 1); 
dissipation3 = NaN(n_profiles,1); 

for ensemble = 1:n_profiles
    Deez(ensemble,:) = AQHR_StructureFunction(vstar(starts(ensemble): ...
        starts(ensemble) + ens_length - 1,:), r_min, r_max);
end

r_struct = 0.022*(r_min:r_max);
r_struct = r_struct.^(2/3); 
fitStruct = NaN(n_profiles,2);

for ensemble = 1:n_profiles
    fitStruct = robustfit(r_struct, Deez(ensemble,:),'ols'); 
    %if fitStruct > 0
        dissipation3(ensemble) = (abs(fitStruct(2)/2.1)).^(3/2);
    %end
end

k3 = sgolayfilt(log10(dissipation3),2,35);

plot(diss_time, k3,'k');

%% GUERRA & THOMSON'S METHOD

dissipation4 = NaN(n_profiles,1);

wavenums = (1:floor(size(vstar,2)/2-1))/(0.022*size(vstar,2));
wavefit = wavenums.^(5/3); wavefit2 = wavefit.^(-1);


for ensemble = 1:n_profiles
    dissipation4(ensemble) = Guerra_Dissipation(vstar(starts(ensemble):...
        starts(ensemble)+ens_length-1,:), 0.022, 30, 0.4, wavenums, wavefit);
end

k4 = sgolayfilt(log10(dissipation4),2,35);
plot(diss_time,k4,'c');

legend('Sreenivasan 1995','Veron & Melville','Wiles et al 2006',...
    'Guerra & Thomson 2011');



