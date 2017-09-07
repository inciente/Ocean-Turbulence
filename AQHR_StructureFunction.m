<<<<<<< HEAD
function Dmat = AQHR_StructureFunction(vstar, r_min, r_max)
=======
function [D] = AQHR_StructureFunction(vstar, r_min, r_max)
>>>>>>> f7e5b6d384b3025f8909a6a2f5aaf5d1b1418ed0

%Vstar is a matrix of along-beam velocity measurements where the vertical 
%dimension is time and the horizontal is bins.
%
%r_min and r_max must be given as separation between bins for which
%structure functions will be computed and average. 
%
%Note that setting a large r_max may cause results to diverge from the r^(2/3) 
%trend, since there will be few velocity differences to average for large r.
%
%vstar = detrend(vstar')';
<<<<<<< HEAD
%vstar = vstar - nanmean(vstar,2);
vstar = vstar - nanmean(vstar, 1);

if r_min == 1
    warning('Wiles et al (2006) suggests two consecutive bins should not be compared', ...
        'Please set r_min > 1')
end
=======
vstar = vstar - nanmean(vstar, 1);
%According to Wiles et al (2006) we shouldn't compare consecutive bins, so 
%r_min is usually set to 2.
r_max = floor(r_max); r_min = floor(r_min); 
>>>>>>> f7e5b6d384b3025f8909a6a2f5aaf5d1b1418ed0


%Matrix form of the structure function, before temporal averaging
Dmat = zeros(size(vstar,1),r_max - r_min + 1);


%First we want to identify which values of r (distances between bins) we
%can test and then how to find the relevant observations in vstar to make
%those comparisons.
for r = r_min:r_max
    
    %How many times can this r be measured inside our profile??
    cases = size(vstar,2) - r;
<<<<<<< HEAD
=======
    %Entries Dmat(time,r) will be the average of the squared differences of
    %along-beam velocity at separations r at time t.
    
>>>>>>> f7e5b6d384b3025f8909a6a2f5aaf5d1b1418ed0
    %Dmat(:,r-r_min+1) = mean((vstar(:,1:cases) - vstar(:,r+1:r+cases)).^2,2);
    for k = 1:cases;
        
        %For all times (indicated by : in Dmat), we store the sum of
        %squared differences between bins separated by r
        Dmat(:,r-r_min+1) = Dmat(:,r-r_min+1) + (vstar(:,k) - vstar(:,k+r)).^2; 
    end
    Dmat(:,r-r_min+1) = Dmat(:,r-r_min+1)/cases;
end

<<<<<<< HEAD
Dmat = nanmean(Dmat,1); %We average over time.
=======
D = nanmean(Dmat,1); %We average over time.
>>>>>>> f7e5b6d384b3025f8909a6a2f5aaf5d1b1418ed0

end
