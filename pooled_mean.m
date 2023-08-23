function [pooled_pdfs,POOLED_AGE] = pooled_mean(period_no,pdfs,P,delta)
%% function for calculating the pooled mean of modeled ages within a period 
%INPUT
%period_no: period number 
%pdfs: vector containing estimated PDF of the ages
%P: structure containing the information of all radiocarbon ages
%delta: nearest years to be rounded
%OUTPUT
%pooled_pdfs: matrix containing pooled mean pdfs and higher probabilities
%POOLED_AGE: structure wrapping up the pooled mean age and hpd regions
%%
M = length(P);  %total number of periods
N = zeros(1,M); %listing the number of ages in each period
for i = 1:M
    N(i) = length(P(i).age); 
end
%% calculating pooled mean pdf of ages in the period
ind_b_theta = cumsum([1 N(1:end-1)+2]) + 1; %index of beginning of theta in the sequence
ind_e_theta = cumsum([1 N(1:end-1)+2]) + N; %index of end of theta in the sequence
cal_age = pdfs(:,1);
x = pdfs(:,2:end);
prob = sum(x(:,ind_b_theta(period_no):ind_e_theta(period_no)),2)/N(period_no);
prob = prob/sum(prob);% normalizing to 1
%%
%% caulcating the 68.2% and 95.4% higher posterior probability and hpd regions.
% MatCal 3.0 (2020-08-13)
% Function for 14C age calibration using Bayesian higher posterior
% density analysis of a probability density function of calibrated age.
%
% Please see manuscript for more detailed information:
% Lougheed, B.C. & Obrochta, S.P. (2016). MatCal: Open Source Bayesian
% 14C Age Calibration in Matlab. Journal of Open Research Software. 4(1),
% p.e42. DOI: http://doi.org/10.5334/jors.130
%%
hp_68p2 = zeros(length(prob),1);
hp_95p4 = zeros(length(prob),1);
calprob = [cal_age prob];
hpd = calprob(:,1:2);
hpd = sortrows(hpd,2);
hpd(:,3) = cumsum(hpd(:,2));
% one-sigma highest posterior probability
hpd68_2 = hpd(hpd(:,3) >= 1-erf(1/sqrt(2)),:);
hpd68_2 = sortrows(hpd68_2,1);
id0 = ismember(cal_age,hpd68_2(:,1)); 
hp_68p2(id0) = hpd68_2(:,2);
% two-sigma highest posterior probability
hpd95_4 = hpd(hpd(:,3) >= 1-erf(2/sqrt(2)),:);
hpd95_4 = sortrows(hpd95_4,1);
id1 = ismember(cal_age,hpd95_4(:,1)); 
hp_95p4(id1) = hpd95_4(:,2);
% one-sigma HPD region
ind1 = find(diff(hpd68_2(:,1)) > delta);
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
% two-sigma HPD region
ind2 = find(diff(hpd95_4(:,1)) > delta);
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
%%
%% calculating age at median probability
[~, median_ind] = min(abs(cumsum(calprob(:,2))-0.5));
median_age = round(calprob(median_ind(1),1));
%% reporting results
POOLED_AGE = struct('P68_2_regions',[],'P95_4_regions',[],'Median_age',[]);
POOLED_AGE.P68_2_regions = p68_2;
POOLED_AGE.P95_4_regions = p95_4;
POOLED_AGE.Median_age = median_age; 
pooled_pdfs = [cal_age prob hp_68p2 hp_95p4];
end