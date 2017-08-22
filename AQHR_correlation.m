function [Rcorr, vmean] = AQHR_correlation(vstar, r_min, r_max, meanremoval)
    
    %Compute the spatial correlation function of doppler current profiles, 
    %and the mean squared velocity of an ensemble to normalize power
    %spectra following Veron & Melville (1999)
    
    %Input fields are:
    %vstar - matrix of along-beam velocities with dimensions (time, bins)
    %r_min - the minimum value of r for which Rcorr shall be computed
    %r_max - maximum value of r for which Rcorr shall be computed
    %meanremoval - boolean to tell whether the temporal mean velocity
    %should be removed from each bin's data.
    
    %Do you want to remove the temporal mean velocity from each bin?
    if meanremoval == 1
        vstar = vstar - nanmean(vstar,1); 
    end
    
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
    
end
