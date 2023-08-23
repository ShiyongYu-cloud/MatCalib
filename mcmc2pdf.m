function [pdfs,hpd68_2,hpd95_4] = mcmc2pdf(mcmcsamples,A,B,delta,year_scale)
%% function for calculating the pdf of calendar ages and save results to a file 
%INPUT
%mcmcsamples: randomly generated calendar ages of the raddiocaron ages and boundaries
%A: early boundary of the age sequence
%B: late boundary of the age sequence
%delta: resolution of the age points in years
%year_scale: scale of year(BCE/CE or BP)
%OUTPUT
%pdfs: vector containing estimated PDF at the age points
%hp_68p2: vector containing the 68.2% highest probability 
%hp_95p4: vector containing the 95.4% highest probability
%% preallocate size
if strcmpi(year_scale,'BCE/CE') == 1
   x = A:B;
else
   x = B:A;
end
%resample at delta-yr intervals
x = x(mod(x,delta) == 0);
M = length(x);
N = size(mcmcsamples,2);
pdfs = zeros(M,N+1);
hpd68_2 = zeros(M,N+1);
hpd95_4 = zeros(M,N+1);
%% get results
for i = 1:N
    [cal_age,prob,hp_68p2,hp_95p4] = epdf(mcmcsamples(:,i),A,B,delta,year_scale);
    pdfs(:,i+1) = prob;
    hpd68_2(:,i+1) = hp_68p2;
    hpd95_4(:,i+1) = hp_95p4;
end
pdfs(:,1) = cal_age;
hpd68_2(:,1) = cal_age; 
hpd95_4(:,1) = cal_age;
return;
%%
function [cal_age,prob,hp_68p2,hp_95p4] = epdf(X,A,B,delta,year_scale)
%% function for estimating the empirical pdf and highest probability from MCMC runs
%INPUT
%X: vector specifying random variables obtained from MCMC run 
%A: early bound of the age sequence
%B: late bound of the age sequence
%delta: resolution in years
%year_scale: scale of year(BCE/CE or BP)
%
%OUTPUT
%cal_age: calendar age points at which PDF was estimated 
%prob: vector specifying estimated PDF at age points
%hp_68p2: vector specifying the 68.2% highest posterior probability
%hp_95p4: vector specifying the 95.4% highest posterior probability
%% building up an empirical pdf
%generate points at which pdf will be estimated
if strcmpi(year_scale,'BCE/CE') == 1
   cal_age = A:B;
else
   cal_age = B:A;
end
%resample at delta-yr interval
cal_age = cal_age(mod(cal_age,delta) == 0);
cal_age = cal_age';
M = length(cal_age);
N = length(X);
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
        if X(j)<=cal_age(i)   % Definition of CDF
            p = p + 1;
        else
            q = q + 1;
        end
    end
    F(i) = p/(p + q);   % Calulating Probability
end
% Estimating PDF by differentiation of CDF
for k = 1:M
    fxph = interp1(cal_age,F,cal_age(k) + h,'spline');  % Interpolating value of F(x+h)
    fxmh = interp1(cal_age,F,cal_age(k) - h,'spline');  % Interpolating value of F(x-h)
    prob(k) = (fxph - fxmh)/(2*h); 
    if prob(k) < 0
       prob(k) = 0;
    end        
end                                         
prob = smooth(prob);  
% Normalizing to 1
prob = prob./sum(prob);
%% calculate higher probabilities
calprob = [cal_age prob];
hpd = calprob(:,1:2);
hpd = sortrows(hpd, 2);
hpd(:,3) = cumsum(hpd(:,2));
%1 sigma prob
hpd68_2 = hpd(hpd(:,3) >= 1-erf(1/sqrt(2)), :);
hpd68_2 = sortrows(hpd68_2,1);
id0 = ismember(cal_age,hpd68_2(:,1)); 
hp_68p2(id0) = hpd68_2(:,2);
%2 sigma prob
hpd95_4 = hpd(hpd(:,3) >= 1-erf(2/sqrt(2)), :);
hpd95_4 = sortrows(hpd95_4,1);
id1 = ismember(cal_age,hpd95_4(:,1)); 
hp_95p4(id1) = hpd95_4(:,2);
return;