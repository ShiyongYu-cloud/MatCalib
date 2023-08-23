function [theta,alpha,beta] = initialization(P,A,B,year_scale,CalCurves,boundary_relationship)
%% function for initializing model parameters
%INPUT
%P: structure containing laboratory data of the periods of a sequence
%A: early boundary of the age sequence
%B: late boundary of the age sequence
%year_scale: scale of year ("BCE/CE" or "BP")
%CalCurves: structure containing the extracted calibration curves in the F14C space
%boundary_relationship: vector containing the relationship of two periods ('overlapping',disjoint','contiguous')
%OUTPUT
%theta: structure containing the calendar ages of the periods
%alpha: vector containing the calendar age of early boundary of all periods
%beta: vector containing the calendar ages of late boundary of all periods
%% obtain calibrated ages of the radiocarbon ages in the periods
M = length(P); % number of periods
alpha = zeros(1,M);
beta = zeros(1,M);
theta = struct('age',[]);
for i = 1:M
    x = calib(P(i).f14C,P(i).ferr,P(i).cal_curve,CalCurves,year_scale);  
    if i == 1 % the first period 
        if strcmpi(year_scale,'BCE/CE') == 1 
            true = all(x > A) && all(diff(x)~=0); % obtain unique ages in a period 
            while(~true)
                 x = calib(P(i).f14C,P(i).ferr,P(i).cal_curve,CalCurves,year_scale);  
                 true = all(x > A) && all(diff(x)~=0);
            end
            if strcmpi(P(i).ordering,'unordered') == 1
                theta(i).age = x; % no ordering               
            elseif strcmpi(P(i).ordering,'coeval') == 1
                theta(i).age = repmat(round(mean(x)),1,length(x)); % make all(diff(theta.age))=0    
            elseif strcmpi(P(i).ordering,'ordered') == 1
                theta(i).age = sort(x,'ascend'); % make all(diff(theta.age))>0
            end    
        elseif strcmpi(year_scale,'BP') == 1
            true = all(x < A) && all(diff(x)~=0);
            while(~true)
                x = calib(P(i).f14C,P(i).ferr,P(i).cal_curve,CalCurves,year_scale);  
                true = all(x < A) && all(diff(x)~=0);
            end
            if strcmpi(P(i).ordering,'unordered') == 1
                theta(i).age = x; % no ordering                 
            elseif strcmpi(P(i).ordering,'coeval') == 1
                theta(i).age = repmat(round(mean(x)),1,length(x)); % make all(diff(theta.age))=0     
            elseif strcmpi(P(i).ordering,'ordered') == 1
                theta(i).age = sort(x,'descend'); % make all(diff(theta.age))<0
            end
        end
    elseif i > 1 && i < M %internal periods    
        if strcmpi(year_scale,'BCE/CE') == 1
            if strcmpi(boundary_relationship(i-1),'overlapping') == 1
                true = min(theta(i-1).age) < min(x) && min(x) < max(theta(i-1).age) && max(theta(i-1).age) < max(x) && all(diff(x)~=0);
            else %contiguous and disjoint periods    
                true = max(theta(i-1).age) < min(x) && all(diff(x)~=0);
            end    
            while(~true)
                x = calib(P(i).f14C,P(i).ferr,P(i).cal_curve,CalCurves,year_scale);  
                if strcmpi(boundary_relationship(i-1),'overlapping') == 1
                    true = min(theta(i-1).age) < min(x) && min(x) < max(theta(i-1).age) && max(theta(i-1).age) < max(x) && all(diff(x)~=0);
                else %contiguous and disjoint periods    
                    true = max(theta(i-1).age) < min(x) && all(diff(x)~=0);
                end    
            end 
            if strcmpi(P(i).ordering,'unordered') == 1
                theta(i).age = x; % no ordering                 
            elseif strcmpi(P(i).ordering,'coeval') == 1
                theta(i).age = repmat(round(mean(x)),1,length(x)); % make all(diff(theta.age))=0
            elseif strcmpi(P(i).ordering,'ordered') == 1  
                theta(i).age = sort(x,'ascend'); % make all(diff(theta.age))>0
            end    
        elseif strcmpi(year_scale,'BP') == 1
            if strcmpi(boundary_relationship(i-1),'overlapping') == 1
                true = max(theta(i-1).age) > max(x) && max(x) > min(theta(i-1).age) && min(theta(i-1).age) > min(x) && all(diff(x)~=0);
            else % dcontiguous and disjoint periods
                true = min(theta(i-1).age) > max(x) && all(diff(x)~=0);
            end
            while(~true)
                x = calib(P(i).f14C,P(i).ferr,P(i).cal_curve,CalCurves,year_scale);  
                if strcmpi(boundary_relationship(i-1),'overlapping') == 1
                    true = max(theta(i-1).age) > max(x) && max(x) > min(theta(i-1).age) && min(theta(i-1).age) > min(x) && all(diff(x)~=0);
                else % dcontiguous and disjoint periods
                    true = min(theta(i-1).age) > max(x) && all(diff(x)~=0);
                end
            end 
            if strcmpi(P(i).ordering,'unordered') == 1
                theta(i).age = x; % no ordering                 
            elseif strcmpi(P(i).ordering,'coeval') == 1
                theta(i).age = repmat(round(mean(x)),1,length(x)); % make all(diff(theta.age))=0
            elseif strcmpi(P(i).ordering,'ordered') == 1
                theta(i).age = sort(x,'descend'); % make all(diff(theta.age))<0
            end    
        end 
    else % the last period
        if strcmpi(year_scale,'BCE/CE') == 1
            if strcmpi(boundary_relationship(i-1),'overlapping') == 1
                true = min(theta(i-1).age) < min(x) && min(x) < max(theta(i-1).age) && max(theta(i-1).age) < max(x) && max(x) < B && all(diff(x)~=0);
            else %contiguous and disjoint periods    
                true = max(theta(i-1).age) < min(x) && max(x) < B && all(diff(x)~=0);
            end    
            while(~true)
                x = calib(P(i).f14C,P(i).ferr,P(i).cal_curve,CalCurves,year_scale);  
                if strcmpi(boundary_relationship(i-1),'overlapping') == 1
                    true = min(theta(i-1).age) < min(x) && min(x) < max(theta(i-1).age) && max(theta(i-1).age) < max(x) && max(x) < B && all(diff(x)~=0);
                else %contiguous and disjoint periods    
                    true = max(theta(i-1).age) < min(x) && max(x) < B && all(diff(x)~=0);
                end    
            end
            if strcmpi(P(i).ordering,'unordered') == 1
                theta(i).age = x; % no ordering                 
            elseif strcmpi(P(i).ordering,'coeval') == 1
                theta(i).age = repmat(round(mean(x)),1,length(x)); % make all(diff(theta.age))=0
            elseif strcmpi(P(i).ordering,'ordered') == 1
                theta(i).age = sort(x,'ascend'); % make all(diff(theta.age))>0
            end    
        elseif strcmpi(year_scale,'BP') == 1
            if strcmpi(boundary_relationship(i-1),'overlapping') == 1
                true = max(theta(i-1).age) > max(x) && max(x) > min(theta(i-1).age) && min(theta(i-1).age) > min(x) && min(x) > B && all(diff(x)~=0);
            else % dcontiguous and disjoint periods
                true = min(theta(i-1).age) > max(x) && min(x) > B && all(diff(x)~=0);    
            end
            while(~true)
                x = calib(P(i).f14C,P(i).ferr,P(i).cal_curve,CalCurves,year_scale);  
                if strcmpi(boundary_relationship(i-1),'overlapping') == 1
                    true = max(theta(i-1).age) > max(x) && max(x) > min(theta(i-1).age) && min(theta(i-1).age) > min(x) && min(x) > B && all(diff(x)~=0);
                else % dcontiguous and disjoint periods
                    true = min(theta(i-1).age) > max(x) && min(x) > B && all(diff(x)~=0);    
                end
            end
            if strcmpi(P(i).ordering,'unordered') == 1
                theta(i).age = x; % no ordering                 
            elseif strcmpi(P(i).ordering,'coeval') == 1
                theta(i).age = repmat(round(mean(x)),1,length(x)); % make all(diff(theta.age))=0
            elseif strcmpi(P(i).ordering,'ordered') == 1
                theta(i).age = sort(x,'descend'); % make all(diff(theta.age))<0
            end    
        end 
    end
end
%% obtain the early boundary of the periods
for i = 1:M
    if i == 1 % the first period
        if strcmpi(year_scale,'BCE/CE') == 1
            x = (A + min(theta(i).age))/2;
            alpha(i) = trirnd(x,min(theta(i).age),min(theta(i).age),1);
        elseif strcmpi(year_scale,'BP') == 1
            x = (A + max(theta(i).age))/2;
            alpha(i) = trirnd(max(theta(i).age),max(theta(i).age),x,1);
        end
    else %other periods
        if strcmpi(year_scale,'BCE/CE') == 1
            if strcmpi(boundary_relationship(i-1),'overlapping') == 1
                id = theta(i-1).age < min(theta(i).age);
                x = (max(theta(i-1).age(id)) + min(theta(i).age))/2;
                alpha(i) = trirnd(max(theta(i-1).age(id)),x,min(theta(i).age),1);
            elseif strcmpi(boundary_relationship(i-1),'contiguous') == 1
                x = (max(theta(i-1).age) + min(theta(i).age))/2;
                alpha(i) = trirnd(max(theta(i-1).age),x,min(theta(i).age),1);
            elseif strcmpi(boundary_relationship(i-1),'disjoint') == 1
                x = (max(theta(i-1).age) + min(theta(i).age))/2;
                x = (x + min(theta(i).age))/2;
                alpha(i) = trirnd(x,min(theta(i).age),min(theta(i).age),1);
            end
        elseif strcmpi(year_scale,'BP') == 1
            if strcmpi(boundary_relationship(i-1),'overlapping') == 1
                id = theta(i-1).age > max(theta(i).age);
                x = (min(theta(i-1).age(id)) + max(theta(i).age))/2;
                alpha(i) = trirnd(max(theta(i).age),x,min(theta(i-1).age(id)),1);
            elseif strcmpi(boundary_relationship(i-1),'contiguous') == 1
                x = (min(theta(i-1).age) + max(theta(i).age))/2;
                alpha(i) = trirnd(max(theta(i).age),x,min(theta(i-1).age),1);
            elseif strcmpi(boundary_relationship(i-1),'disjoint') == 1
                x = (min(theta(i-1).age) + max(theta(i).age))/2;
                x = (x + max(theta(i).age))/2;
                alpha(i) = trirnd(max(theta(i).age),max(theta(i).age),x,1);
            end
        end 
    end
    alpha(i) = round(alpha(i));
end   
%% obtain the late boundary of the periods
for i = 1:M
    if i < M % first M-1 periods
        if strcmpi(year_scale,'BCE/CE') == 1
            if strcmpi(boundary_relationship(i),'overlapping') == 1
                id = theta(i+1).age > max(theta(i).age);
                x = (max(theta(i).age) + min(theta(i+1).age(id)))/2;
                beta(i) = trirnd(max(theta(i).age),x,min(theta(i+1).age(id)),1);
            elseif strcmpi(boundary_relationship(i),'contiguous') == 1
                beta(i) = alpha(i+1);
            elseif strcmpi(boundary_relationship(i),'disjoint') == 1
                x = (max(theta(i).age) + min(theta(i+1).age))/2;
                x = (max(theta(i).age) + x)/2;
                beta(i) = trirnd(max(theta(i).age),max(theta(i).age),x,1);
            end
        elseif strcmpi(year_scale,'BP') == 1
            if strcmpi(boundary_relationship(i),'overlapping') == 1
                id = theta(i+1).age < min(theta(i).age);
                x = (min(theta(i).age) + max(theta(i+1).age(id)))/2;
                beta(i) = trirnd(max(theta(i+1).age(id)),x,min(theta(i).age),1);
            elseif strcmpi(boundary_relationship(i),'contiguous') == 1
                beta(i) = alpha(i+1);
            elseif strcmpi(boundary_relationship(i),'disjoint') == 1
                x = (min(theta(i).age) + max(theta(i+1).age))/2;
                x = (min(theta(i).age) + x)/2;
                beta(i) = trirnd(x,min(theta(i).age),min(theta(i).age),1);
            end
        end 
    else % the last period
        if strcmpi(year_scale,'BCE/CE') == 1
            x = (B + max(theta(i).age))/2;
            beta(i) = trirnd(max(theta(i).age),max(theta(i).age),x,1);
        elseif strcmpi(year_scale,'BP') == 1
            x = (B + min(theta(i).age))/2;
            beta(i) = trirnd(x,min(theta(i).age),min(theta(i).age),1);
        end
    end
    beta(i) = round(beta(i));
end  
return;
%%
function cal_age = calib(F_14C_val,F_14C_err,cal_curves,CalCurves,year_scale)
%% Function for calibrating 14C ages in a period using Bayesian highest posterior density analysis
%INPUT
% F_14C: laboratory 14C age in F14C space
% F_err: laboratory 14C age error in F14C space
% cal_curves: calibration curves associated with the 14C ages 
% CalCurve: structure containing the calibration curves
% year_scale: scale of year to be reported (e.g. BC/AD or BP)
%OUTPUT
% cal_age: vector containing the calibrated ages
%%
N = length(F_14C_val);
cal_age = zeros(1,N);
interpres = 1; % interpolate the curve to annual resolution
for i = 1:N
    % retrieve the calibration dataset
    if strcmpi(cal_curves(i),'intcal13') == 1
        curve_cal_age = CalCurves.intcal13(:,1);
        curve_F14_val = CalCurves.intcal13(:,2);
        curve_F14_err = CalCurves.intcal13(:,3);
        hi_curve_cal_age = curve_cal_age(1):interpres:curve_cal_age(end);
        hi_curve_F14_val = interp1(curve_cal_age, curve_F14_val, hi_curve_cal_age);
        hi_curve_F14_err = interp1(curve_cal_age, curve_F14_err, hi_curve_cal_age);
    elseif strcmpi(cal_curves(i),'marine13') == 1
        curve_cal_age = CalCurves.marine13(:,1);
        curve_F14_val = CalCurves.marine13(:,2);
        curve_F14_err = CalCurves.marine13(:,3);
        hi_curve_cal_age = curve_cal_age(1):interpres:curve_cal_age(end);
        hi_curve_F14_val = interp1(curve_cal_age, curve_F14_val, hi_curve_cal_age);
        hi_curve_F14_err = interp1(curve_cal_age, curve_F14_err, hi_curve_cal_age);
    elseif strcmpi(cal_curves(i),'shcal13') == 1
        curve_cal_age = CalCurves.shcal13(:,1);
        curve_F14_val = CalCurves.shcal13(:,2);
        curve_F14_err = CalCurves.shcal13(:,3);   
        hi_curve_cal_age = curve_cal_age(1):interpres:curve_cal_age(end);
        hi_curve_F14_val = interp1(curve_cal_age, curve_F14_val, hi_curve_cal_age);
        hi_curve_F14_err = interp1(curve_cal_age, curve_F14_err, hi_curve_cal_age);
    elseif strcmpi(cal_curves(i),'intcal20') == 1
        curve_cal_age = CalCurves.intcal20(:,1);
        curve_F14_val = CalCurves.intcal20(:,2);
        curve_F14_err = CalCurves.intcal20(:,3);
        hi_curve_cal_age = curve_cal_age(1):interpres:curve_cal_age(end);
        hi_curve_F14_val = interp1(curve_cal_age, curve_F14_val, hi_curve_cal_age);
        hi_curve_F14_err = interp1(curve_cal_age, curve_F14_err, hi_curve_cal_age);
    elseif strcmpi(cal_curves(i),'marine20') == 1
        curve_cal_age = CalCurves.marine20(:,1);
        curve_F14_val = CalCurves.marine20(:,2);
        curve_F14_err = CalCurves.marine20(:,3);
        hi_curve_cal_age = curve_cal_age(1):interpres:curve_cal_age(end);
        hi_curve_F14_val = interp1(curve_cal_age, curve_F14_val, hi_curve_cal_age);
        hi_curve_F14_err = interp1(curve_cal_age, curve_F14_err, hi_curve_cal_age);
    elseif strcmpi(cal_curves(i),'shcal20') == 1
        curve_cal_age = CalCurves.shcal20(:,1);
        curve_F14_val = CalCurves.shcal20(:,2);
        curve_F14_err = CalCurves.shcal20(:,3);        
        hi_curve_cal_age = curve_cal_age(1):interpres:curve_cal_age(end);
        hi_curve_F14_val = interp1(curve_cal_age, curve_F14_val, hi_curve_cal_age);
        hi_curve_F14_err = interp1(curve_cal_age, curve_F14_err, hi_curve_cal_age);
    end
    % calculate the highest posterior probabilities
    a = (F_14C_val(i) - hi_curve_F14_val).^2;
    b = 2*(F_14C_err(i).^2 + hi_curve_F14_err.^2);
    c = sqrt(F_14C_err(i).^2 + hi_curve_F14_err.^2);
    cal_prob = exp(-a./b)./c;
    % normalize the probabilities to 1
    cal_prob = cal_prob/sum(cal_prob);
    % draw a calendar age from the pdf
    cal_age(i) = durand(hi_curve_cal_age,cal_prob,1,1);
end
%% convert to BC/AD scale if needed
if strcmpi(year_scale,'BCE/CE') == 1
   cal_age = 1950 - cal_age;
end
return;
%%
function x = durand(values,prob,M,N)
%% Function for drawing a random number from a list at given probability
values = values(:);
prob = prob(:);
if sum(prob)~=1
   prob = prob/sum(prob);
end
L = length(prob);
K = M*N;
psup = cumsum(prob);
pinf = [0; psup(1:end-1)];
Pinf = kron(ones(1,K),pinf(:));
Psup = kron(ones(1,K),psup(:));
u = rand(1,K);
U = kron(ones(L,1),u);
C = (U>Pinf) & (U<Psup);
V = kron(values(:),ones(1,K));
X = V.*C;
x = sum(X);
x = reshape(x,M,N);
return;