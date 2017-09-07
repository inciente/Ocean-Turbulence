<<<<<<< HEAD
function [spectra] = AQHR_correlation(vstar, r_min, r_max, ...
    meanremoval, win_length, overlap, bindistance)
=======
function [Rcorr, vmean] = AQHR_correlation(vstar, r_min, r_max, meanremoval)
>>>>>>> f7e5b6d384b3025f8909a6a2f5aaf5d1b1418ed0
    
    %Compute the spatial correlation function of doppler current profiles, 
    %and the mean squared velocity of an ensemble to normalize power
    %spectra following Veron & Melville (1999)
    
<<<<<<< HEAD
    %Then find the average power density spectrum from all profiles in the
    %ensemble.
    
=======
>>>>>>> f7e5b6d384b3025f8909a6a2f5aaf5d1b1418ed0
    %Input fields are:
    %vstar - matrix of along-beam velocities with dimensions (time, bins)
    %r_min - the minimum value of r for which Rcorr shall be computed
    %r_max - maximum value of r for which Rcorr shall be computed
    %meanremoval - boolean to tell whether the temporal mean velocity
    %should be removed from each bin's data.
<<<<<<< HEAD
    %win_length - window length to apply the Welch spectral method
    %overlap - overlap of windows in Welch's method
=======
>>>>>>> f7e5b6d384b3025f8909a6a2f5aaf5d1b1418ed0
    
    %Do you want to remove the temporal mean velocity from each bin?
    %McLane-mounted aquadopps are bound to produce high self correlation if
    %mean velocities are not substracted from each temporal profile, since
    %these velocities may dominate the measurements and make them all
    %fairly similar.
    if meanremoval == 1
        vstar = vstar - nanmean(vstar,2); 
    end
    
<<<<<<< HEAD
    %vmean = mean(vstar.^2,2); 
    
    %Create an array to store the correlation function
    %Rcorr = zeros(size(vstar,1),r_max - r_min +1); 
    %spectra = zeros(size(vstar,1), floor((r_max - r_min +1)/2)-1); 
    spectra = zeros(size(vstar,1), floor(win_length/2));
    
    %Compute the value of the correlation function at all distances r
%     for r = r_min:r_max
%         
%         %How many measurements of this r's correlation will we get to
%         %average?
%         cases = size(vstar,2) - r;
%         
%         for k=1:cases
%             Rcorr(:,r - r_min + 1) = Rcorr(:,r-r_min+1) + vstar(:,k).*vstar(:,k+r);
%         end
%         Rcorr(:,r - r_min + 1) = Rcorr(:,r-r_min+1)/cases;
%         %With this, correlation functions have been computed at all times.
%     end
   
    
    for t = 1:size(vstar,1)
    %Compute psd for each correlation function in time
        %one_psd = obmPSpec(vstar(t,r_min:r_max), bindistance, win_length, overlap);
        one_psd = abs(fft(vstar(t,r_min:r_max))).^2;
        %spectra(t,:) = one_psd.psd;
        spectra(t,:) = one_psd(2:size(spectra,2)+1);
    
    end
    
    spectra = spectra./ceil((r_max - r_min));
    spectra = mean(spectra,1);
=======
    vmean = mean(mean(vstar.^2));
    
    %Create an array to store the correlation function
    Rcorr = zeros(size(vstar,1),r_max - r_min +1); 
    
    for r = r_min:r_max
        
        %How many measurements of this r's correlation will we get to
        %average?
        cases = size(vstar,2) - r;
        
        for k=1:cases
            Rcorr(:,r - r_min + 1) = Rcorr(:,r-r_min+1) + vstar(:,k).*vstar(:,k+r);
        end
            Rcorr(:,r - r_min + 1) = Rcorr(:,r-r_min+1)/cases;
    end
    Rcorr = mean(Rcorr,1); 
>>>>>>> f7e5b6d384b3025f8909a6a2f5aaf5d1b1418ed0
    
end
