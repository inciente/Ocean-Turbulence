function [D] = AQHR_StructureFunction(vstar, binheight)

%Vstar is a matrix of along-beam velocity measurements where the vertical 
%dimension is time and the horizontal is bins
%
%Dprofs is the matrix where structure function results are going to be stored.
%
%binheight is used to create the values of r^(2/3) when fitting D 
%
%D, fit and stats are the structure function for the set in question, the 
%values of A and N (see Wiles et al. 2006) and the goodness-of-fit
%statistics respectively.

%Time average of differences in fluctuations u' between the first and all
%other bins:

vstar = vstar - nanmean(vstar); 

D = nanmean((vstar(:,1) - vstar(:,4:end)).^2);

%We now fit 

end



