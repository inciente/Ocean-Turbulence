function spectrum = norm_spectrum(vstar, meanremoval);
    %Compute a normalized spectrum for data in vstar according to the
    %method described by Sreenivasan (1995). 
    
    if meanremoval ==1
       vstar = vstar - nanmean(vstar,2);  
    end
    
    vmean = mean(vstar.^2,2); 
    
    spec_matrix = NaN(size(vstar,1), floor(size(vstar,2)/2-1)); 
    
    for k=1:size(vstar,1)
        spectrum = abs(fft(vstar(k,:)));
        spectrum = spectrum(2:size(spec_matrix,2)+1).^2;
        spec_matrix(k,:) = spectrum;
        
    end
    
    %spec_matrix = spec_matrix./(size(vstar,2)-1);
    sums = sum(spec_matrix,2); 
    spec_matrix = vmean.*spec_matrix./sums;
    
    spectrum = mean(spec_matrix,1); 
end
