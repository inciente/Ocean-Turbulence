clc; clear all

%Script to visualize statistics of results looking for the wavenumber
%bounds of the inertial subrange in AQDP measurements. 
%Since the end goal of this entire project is to generate good estimates of
%epsilon, sensibility to range is also evaluated against Guerra & Thomson's
%methodology, where they cherrypick estimates.

%Load results from the sixth section of Master_AQHR_epsilon.m
load('k_range.mat');

%Let's start by generating a boxplot matrix of RMSError between spectra and
%fits for various wavenumber ranges

high_bound = NaN(705*9,length(hi_end));

figure;

abs_start = 1;
for element = 1:length(lo_end)
    %Element fixes the lower bound of the wavenumber range, and for each
    %we're going to plot all tested higher bounds (only if hi > lo+1)
    start = find(hi_end == lo_end(element)+2);
    range = abs_start + (length(hi_end) - start);
    
    subplot(3,3,element)
    boxplot(RMSE_data(:,abs_start:range), hi_end(start:end));
    title(['min(k) = ',num2str(lo_end(element))])
    ylim([0 0.5e-6])
    grid on
    
    high_bound((element-1)*705+1:element*705,start:end)... 
        = RMSE_data(:,abs_start:range);
    
    
    abs_start = range+1; 
end



