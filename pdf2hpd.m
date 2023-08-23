function [ALPHA,BETA,THETA,CAGES] = pdf2hpd(pdfs,delta,P,year_scale) 
%% function for calculating the highest posterior density (HPD) regions
%INPUT 
%pdfs: estimated pdf of parameters
%delta: nearest years to be rounded
%P: structure containing the information of all radiocarbon ages
%year_scale: scale of year (BCE/CE or BP)
%OUTPUT
%ALPHA: structure containing the 68.2%, %95.4% highest posterior density regions, and the median prabability age 
%BETA:  structure containing the 68.2%, %95.4% highest posterior density regions, and the median prabability age 
%THETA: structure containing the 68.2%, %95.4% highest posterior density regions, and the median prabability age 
%CAGES: structure wrapping up ALPHA, BETA, and THETA
%%
M = length(P);  %total number of periods
N = zeros(1,M); %listing the number of ages in each period
for i = 1:M
    N(i) = length(P(i).age); 
end
ALPHA = struct('Period_number',[],'Boundary',[],'P68_2_regions',[],'P95_4_regions',[],'Median_prob_age',[],'Year_type',[]);
BETA  = struct('Period_number',[],'Boundary',[],'P68_2_regions',[],'P95_4_regions',[],'Median_prob_age',[],'Year_type',[]);
THETA = struct('Period_number',[],'Lab_code',[],'Radiocarbon_age',[],'P68_2_regions',[],'P95_4_regions',[],'Median_prob_age',[],'Year_type',[]);
CAGES = struct('Period_number',[],'SampleID',[],'P68_2_regions',[],'P95_4_regions',[],'Median_prob_age',[],'Year_type',[]);
cal_age = pdfs(:,1);
prob = pdfs(:,2:end);
ind_alpha = cumsum([1 N(1:end-1)+2]); %index of alpha in the age sequence
ind_b_theta = cumsum([1 N(1:end-1)+2]) + 1; %index of beginning of theta in the sequence 
%ind_e_theta = cumsum([1 N(1:end-1)+2]) + N; %index of end of theta in the sequence
ind_beta = cumsum([1 N(1:end-1)+2]) + N + 1;  %index of beta in the age sequence
ID = cumsum([0 N]);
for i = 1:M
    prob_alpha = prob(:,ind_alpha(i));          %alphas
    ALPHA(i).Period_number = ['P', ' ', num2str(i)];
    ALPHA(i).Boundary = 'early';
    [ALPHA(i).P68_2_regions,ALPHA(i).P95_4_regions,ALPHA(i).Median_prob_age] = hpd(cal_age,prob_alpha,delta);
    ALPHA(i).Year_type = year_scale;
    %
    prob_beta = prob(:,ind_beta(i));            %betas    
    BETA(i).Period_number = ['P', ' ', num2str(i)];
    BETA(i).Boundary = 'late';
    [BETA(i).P68_2_regions,BETA(i).P95_4_regions,BETA(i).Median_prob_age] = hpd(cal_age,prob_beta,delta);
    BETA(i).Year_type = year_scale; 
    %
    CAGES(ind_alpha(i)).Period_number = ['P',' ', num2str(i)];
    CAGES(ind_alpha(i)).SampleID = ['Period',' ', num2str(i),' early boundary'];
    [CAGES(ind_alpha(i)).P68_2_regions,CAGES(ind_alpha(i)).P95_4_regions,CAGES(ind_alpha(i)).Median_prob_age] = hpd(cal_age,prob_alpha,delta);
    CAGES(ind_alpha(i)).Year_type = year_scale;
    CAGES(ind_beta(i)).Period_number = ['P', ' ', num2str(i)];
    CAGES(ind_beta(i)).SampleID = ['Period', ' ', num2str(i),' late boundary'];
    [CAGES(ind_beta(i)).P68_2_regions,CAGES(ind_beta(i)).P95_4_regions,CAGES(ind_beta(i)).Median_prob_age] = hpd(cal_age,prob_beta,delta);
    CAGES(ind_beta(i)).Year_type = year_scale;
    for j = 1:N(i)
        prob_theta = prob(:,ind_b_theta(i)+j-1); %thetas
        %[THETA(i,j).P68_2_HPD_regions,THETA(i,j).P95_4_HPD_regions,THETA(i,j).Median_prob_age] = ci(cal_age,prob_theta,delta);
        %THETA(i,j).Phase_number = ['P' num2str(i)];
        %THETA(i,j).Lab_code = R(i).sampleID{j};
        %THETA(i,j).Year_type = year_scale;
        THETA(ID(i)+j).Period_number = ['P', ' ', num2str(i)];
        THETA(ID(i)+j).Lab_code = P(i).lab_code{j};
        THETA(ID(i)+j).Radiocarbon_age = [num2str(P(i).age(j)),char(177),num2str(P(i).err(j))];
        [THETA(ID(i)+j).P68_2_regions,THETA(ID(i)+j).P95_4_regions,THETA(ID(i)+j).Median_prob_age] = hpd(cal_age,prob_theta,delta);
        THETA(ID(i)+j).Year_type = year_scale;
        %
        CAGES(ind_b_theta(i)+j-1).Period_number = ['P', ' ', num2str(i)];
        CAGES(ind_b_theta(i)+j-1).SampleID = strcat(P(i).lab_code{j},' (',num2str(P(i).age(j)),char(177),num2str(P(i).err(j)),')');
        [CAGES(ind_b_theta(i)+j-1).P68_2_regions,CAGES(ind_b_theta(i)+j-1).P95_4_regions,CAGES(ind_b_theta(i)+j-1).Median_prob_age] = hpd(cal_age,prob_theta,delta);
        CAGES(ind_b_theta(i)+j-1).Year_type = year_scale;
    end
end
return;
%%
function [p68_2,p95_4,median_age] = hpd(cal_age,prob,delta)
%% function for estimating the highest posterior density (hpd) regions of a calendar age 
%INPUT
%cal_age: modeled calendar age points
%prob: probabilities at the calendar age points   
%delta: nearest age to be rounded in years
%OUTPUT
%p68_2: 68.2% highest posterior density regions of the modeled calendar ages 
%p95_4: 95.4% highest posterior density regions of the modeled calendar ages
%median_age: modeld calendar at the median probability
%%
%% Calculating the highest posterior density (HPD) regions
% MatCal 3.0 (2020-08-13)
% Function for 14C age calibration using Bayesian higher posterior
% density analysis of a probability density function of calibrated age.
%
% Please see manuscript for more detailed information:
% Lougheed, B.C. & Obrochta, S.P. (2016). MatCal: Open Source Bayesian
% 14C Age Calibration in Matlab. Journal of Open Research Software. 4(1),
% p.e42. DOI: http://doi.org/10.5334/jors.130
%%
cal_age = cal_age(:);
prob = prob(:);
calprob = [cal_age prob];
hpd = calprob(:,1:2);
hpd = sortrows(hpd,2);
hpd(:,3) = cumsum(hpd(:,2));
% one-sigma highest posterior density regions
hpd68_2 = hpd(hpd(:,3) >= 1-erf(1/sqrt(2)),:);
hpd68_2 = sortrows(hpd68_2,1);
ind1 = find(diff(hpd68_2(:,1)) > delta);
if isempty(ind1) == 1
	p68_2(1,1) = hpd68_2(end,1);
	p68_2(1,2) = hpd68_2(1,1);
	%p68_2(1,3) = sum(hpd68_2(1:end,2));
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
% two-sigma highest posterior density regions
hpd95_4 = hpd(hpd(:,3) >= 1-erf(2/sqrt(2)),:);
hpd95_4 = sortrows(hpd95_4,1);
ind2 = find(diff(hpd95_4(:,1)) > delta);
if isempty(ind2) == 1
	p95_4(1,1) = hpd95_4(end,1);
	p95_4(1,2) = hpd95_4(1,1);
	%p95_4(1,3) = sum(hpd95_4(1:end,2));
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
%% calculate age at median probability
[~, median_ind] = min(abs(cumsum(calprob(:,2))-0.5));
median_age = round(calprob(median_ind(1),1));
return;