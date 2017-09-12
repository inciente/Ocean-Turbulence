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

legend('Sreenivasan 1995','Veron & Melville 1999','Wiles et al 2006',...
    'Guerra & Thomson 2017');

%% WAVENUMBER RANGE TESTING: SPECTRA

wavenums = (1:floor(size(vstar,2)/2))/(0.022*size(vstar,2));
wavefit = wavenums.^(5/3);
% Here we obtain spectra for a few McLane profiles (check start_dt, end_dt)
% and test which is the best wavenumber range to use when fitting ^(-5/3) 
% curves to data to get the least possible amount of noise.
vstar = vstar(starts(1):starts(end)+ens_length-1,:);

%vstar_TKE = vstar - mean(vstar,2); 
%vstar_TKE = vstar - mean(vstar,1); 
vstar_TKE = NaN(size(vstar));

for k = 1:length(starts)
    vstar_TKE((k-1)*ens_length+1:k*ens_length,:) = vstar((k-1)*ens_length+1:...
        k*ens_length,:) - mean(vstar((k-1)*ens_length+1:k*ens_length,:));
end



%Begin by defining the range to be surveyed. 
lo_end = 1:9;
hi_end = 3:15;

%Obtain spectra at all times before making any comparisons
spectra = NaN(n_profiles*ens_length, 16); 
win_length = size(vstar,2);

for time=1:length(spectra)
    
    %Save spectra for a given moment in our new matrix
    spec1 = obmPSpec(vstar_TKE(time,:), 0.022, win_length,0.4);
    spec1 = spec1.psd;
    %Compensate slope in spectra to make them flat
    spectra(time,:) = spec1.*wavefit';
end

%Now get the epsilon estimate that would come out for each ensemble without
%really restricting the inertial subrange wavenumbers.

ref_epsilon = NaN(length(starts), 1); 
vmean = NaN(length(starts),1); 

for ensemble = 1:length(starts)
    
    %Find min(abs(slope)) to compute dissipation from corresponding spectrum
    slope_list = NaN(ens_length,1); 
    vmean(ensemble) = abs(mean(mean(vstar((ensemble-1)*ens_length+1:...
        ensemble*ens_length,:))));
    
    for t_in = 1:ens_length
      
        index = (ensemble-1)*ens_length + t_in;
        fit = robustfit(wavenums, spectra(index,:), 'ols');
        slope_list(t_in) = abs(fit(2)); 
        
    end
    
    %Find min(abs(slope)) and get epsilon as written in Guerra (2017)
    slopeind = find(slope_list == min(slope_list)); 
    slopeind = (ensemble-1)*ens_length + slopeind;
    ref_epsilon(ensemble) = (mean(spectra(slopeind,:))*vmean(ensemble)/...
        (2*pi*0.69))^1.5;
end

clear slopeind slope_list fit index time ensemble t_in spec1

%RMS error of each spectrum and its fit will be saved here. The horizontal
%dimension of this matrix represents all the wavenumber pairs being tested.
%Dimensions will be added as needed with if statement in for loop.
RMSE_data = NaN(length(spectra),1); 
diss_data = NaN(length(spectra),1);
fits = NaN(length(spectra),2);
curr_col = 1;


%Only test when hi_end > lo_end:

for lo=1:length(lo_end)
    lo
    for hi=1:length(hi_end)
        hi
        if hi_end(hi) > lo_end(lo) + 1
            test_k = wavenums(lo_end(lo):hi_end(hi));
            for time = 1:length(spectra)
                fit = robustfit(test_k, spectra(time,lo_end(lo): ...
                    hi_end(hi)), 'ols');
                fits(time,2*curr_col-1:2*curr_col) = fit;
                RMSE_data(time, curr_col) = mean(((fit(1) + fit(2)*test_k) - ...
                    spectra(time,lo_end(lo): hi_end(hi))).^2);
                diss_data(time, curr_col) = mean(spectra(time,lo_end(lo):hi_end(hi))); 
            end
            curr_col = curr_col+1;
        else
            continue
        end
    end
end

alt_average3 = struct('diss_data',diss_data, 'RMSE_data',RMSE_data, 'fits',fits,'spectra',spectra,'vmean',vmean,...
    'ref_epsilon',ref_epsilon');
%save('k_range.mat','diss_data', 'RMSE_data', 'fits','spectra','vmean',...
%    'wavefit','wavenums','hi_end','lo_end','ref_epsilon','fits')
save('k_range3.mat','alt_average3');


%The End
