function [diff_pdfs,AGE_DIFF] = age_difference(x1,x2,bin_size)
%% function for estimating empirical pdf and the highest posterior denstiy (hpd) regions
%% of the difference between two modeled calendar ages
%INPUT
%x1: vector containing the first modeled calendar age points
%x2: vector containing the second modeled calendar age points
%bin_size: size of bins in years for calculating probability of the age difference
%OUTPUT
%diff_pdfs: matrix containing the estimated probability, 68.2%, and 95.4% 
% highest posterior density at the points of age difference
%AGE_DIFF: structure containing the 68.2% and 95.4% hpd regions as well as
%the median probability of age difference 
%% calculating the difference between two modeled calendar ages
x = abs(x1-x2);
a = round(min(x)-5*std(x)); %lower bound of the age differences
if a < 0
    a = 0;
end    
b = round(max(x)+5*std(x)); %upper bound of the age differences
%% building up the empirical pdf
% generate points at which pdf will be estimated
X = a:b;
% resample at delta-yr interval
X = X(mod(X,bin_size) == 0);
X = X';
M = length(X);
N = length(x);
prob = zeros(M,1);
hp_68p2 = zeros(M,1);
hp_95p4 = zeros(M,1);
F = zeros(M,1);
h = 0.000000001;        % Small value closer to zero for evaluating
                        % numerical differentiation.
% Estimating CDF by its definition
for i = 1:M
    p = 0;              % True Probability
    q = 0;              % False Probability
    for j = 1:N
        if x(j) <= X(i)   % Definition of CDF
            p = p + 1;
        else
            q = q + 1;
        end
    end
    F(i) = p/(p + q);   % Calulating Probability
end
% Estimating PDF by differentiating the CDF
for k = 1:M
    fxph = interp1(X,F,X(k) + h,'spline');  % Interpolating value of F(x+h)
    fxmh = interp1(X,F,X(k) - h,'spline');  % Interpolating value of F(x-h)
    prob(k) = (fxph - fxmh)/(2*h); 
    if prob(k) < 0
       prob(k) = 0;
    end        
end                                         
prob = smooth(prob);  
% Normalizing to 1
prob = prob./sum(prob);
%%
%% calculate highest posterior probabilities and hpd regions
% MatCal 3.0 (2020-08-13)
% Function for 14C age calibration using Bayesian higher posterior
% density analysis of a probability density function of calibrated age.
%
% Please see manuscript for more detailed information:
% Lougheed, B.C. & Obrochta, S.P. (2016). MatCal: Open Source Bayesian
% 14C Age Calibration in Matlab. Journal of Open Research Software. 4(1),
% p.e42. DOI: http://doi.org/10.5334/jors.130
%%
calprob = [X prob];
hpd = calprob(:,1:2);
hpd = sortrows(hpd,2);
hpd(:,3) = cumsum(hpd(:,2));
% one-sigma highest posterior probability
hpd68_2 = hpd(hpd(:,3) >= 1-erf(1/sqrt(2)),:);
hpd68_2 = sortrows(hpd68_2,1);
id0 = ismember(X,hpd68_2(:,1)); 
hp_68p2(id0) = hpd68_2(:,2);
% two-sigma highest posterior probability
hpd95_4 = hpd(hpd(:,3) >= 1-erf(2/sqrt(2)),:);
hpd95_4 = sortrows(hpd95_4,1);
id1 = ismember(X,hpd95_4(:,1)); 
hp_95p4(id1) = hpd95_4(:,2);
% one-sigma HPD region
ind1 = find(diff(hpd68_2(:,1)) > bin_size);
if isempty(ind1) == 1
	p68_2(1,1) = hpd68_2(end,1);
	p68_2(1,2) = hpd68_2(1,1);
    p68_2(1,3) = 1; %fraction of enclosed area in the pdf
else
	z = 0;
	for i = 1:length(ind1)
		z = z + 1;
		indy1(z,1) = ind1(i);
		z = z + 1;
		indy1(z,1) = ind1(i)+1;
	end
	indy1 = [ 1 ; indy1; length(hpd68_2(:,1)) ];
	z=0;
	for i = 2:2:length(indy1)
		z = z+1;
		p68_2(z,1) = hpd68_2(indy1(i),1);
		p68_2(z,2) = hpd68_2(indy1(i-1),1);
		p68_2(z,3) = sum(hpd68_2(indy1(i-1):indy1(i),2));
	end
	p68_2 = flipud(p68_2);
    p68_2(:,3) = p68_2(:,3)/sum(p68_2(:,3)); %fraction of each enclosed area in the pdf
end
% keep only the largest region
[~,idd1] = max(p68_2(:,3));
p68_2 = p68_2(idd1,1:2);
p68_2 = fliplr(p68_2);
% two-sigma HPD region
ind2 = find(diff(hpd95_4(:,1)) > bin_size);
if isempty(ind2) == 1
	p95_4(1,1) = hpd95_4(end,1);
	p95_4(1,2) = hpd95_4(1,1);
    p95_4(1,3) = 1; %fraction of enclosed area in the pdf
else
	z = 0;
	for i = 1:length(ind2)
		z = z + 1;
		indy2(z,1) = ind2(i);
		z = z + 1;
		indy2(z,1) = ind2(i)+1;
	end
	indy2 = [ 1 ; indy2; length(hpd95_4(:,1)) ];
	z=0;
	for i = 2:2:length(indy2)
		z = z+1;
		p95_4(z,1) = hpd95_4(indy2(i),1);
		p95_4(z,2) = hpd95_4(indy2(i-1),1);
		p95_4(z,3) = sum(hpd95_4(indy2(i-1):indy2(i),2));
	end
	p95_4 = flipud(p95_4);
    p95_4(:,3) = p95_4(:,3)/sum(p95_4(:,3)); %fraction of enclosed area in the pdf
end
% keep only the largest region
[~,idd2] = max(p95_4(:,3));
p95_4 = p95_4(idd2,1:2);
p95_4 = fliplr(p95_4);
%%
%% calculating age difference at median probability
[~, median_ind] = min(abs(cumsum(calprob(:,2))-0.5));
median_prob_diff = round(calprob(median_ind(1),1));
%% reporting results
AGE_DIFF = struct('P68_2_region',[],'P95_4_region',[],'Median_age_diff',[]);
AGE_DIFF.P68_2_region = p68_2;
AGE_DIFF.P95_4_region = p95_4;
AGE_DIFF.Median_age_diff = median_prob_diff; 
diff_pdfs = [X prob hp_68p2 hp_95p4];
end