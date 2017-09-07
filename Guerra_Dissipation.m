function epsilon = Guerra_Dissipation(vstar, bindistance, win_len, overlap, ...
    wavenums, wavecomp)

    %Computes dissipation through the TKE spectrum of velocity data. This is
    %most similar to the technique used by Guerra & Thomson
    %Wavecomp is an array of wavenumbers k^(5/3) corresponding to the
    %desired spectra, while wavenums is just k.
    vmean = mean(mean(vstar)); 
    %Extract turbulent fluctuations only
    vstar = vstar - mean(vstar,2);
    %Where to store spectra for all profiles in the ensemble vstar
    Spectra = NaN(size(vstar,1), floor(win_len/2)); 
    Slopes = NaN(size(vstar,1),1);
    
    %Wavenumbers encompassing the subinertial range
    k1 = 1;
    k2 = 14;
    
    for time=1:size(vstar,1)
        Spec = obmPSpec(vstar(time,:), bindistance, win_len, overlap);
        Spectra(time,:) = Spec.psd;
        Spectra(time,:) = Spectra(time,:).*wavecomp;
        
        %Now compute the slope of compensated spectra. We want these to be
        %closest to zero as possible.
        slop = robustfit(wavenums, Spectra(time,:),'ols');
        Slopes(time) = abs(slop(2));
    end
    
    %Find which of all fits had the minimum absolute slope
    slop = find(Slopes == min(Slopes));
    epsilon = (mean(Spectra(slop,k1:k2))*abs(vmean)/(2*pi*0.69))^1.5;
    
end
